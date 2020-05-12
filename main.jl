# add Serialization, DataStructures, CoordinateTransformations, Rotations, DataFrames, Missings, Distributions, AngleBetweenVectors, LinearAlgebra, FileIO, Colors
# add https://github.com/yakir12/DungAnalyse.jl
# add AbstractPlotting#use-vertical-dims-from-font CairoMakie#jk/scale_text MakieLayout#master

# 8 -> 1:58
# 4 -> 1:52
# 2 -> 2:49

using DungAnalyse, Serialization, DataStructures, CoordinateTransformations, Rotations, DataFrames, Missings, Distributions, AngleBetweenVectors, LinearAlgebra, StatsBase, OnlineStats
using CairoMakie, MakieLayout, FileIO, AbstractPlotting, Colors
import AbstractPlotting:px
CairoMakie.activate!()
# include("plot_utils.jl")

using PrettyTables, Measurements, HypothesisTests, GLM, MixedModels, DelimitedFiles, Printf

############# data preparation ###################

keep = (:displace_location, :displace_direction , :nest_coverage, :transfer)
discared = ("person", "pellet")
function parsetitle(title, r)
    run = r.data
    d = Dict{Symbol, Any}(k => missing for k in keep)
    d[:nest_coverage] = "open"
    d[:nest] = run.originalnest
    d[:feeder] = run.feeder
    d[:fictive_nest] = run.nest
    d[:track] = run.track
    d[:title] = title
    d[:comment] = r.metadata.comment
    for kv in split(title, ' ')
        k, v = split(kv, '#')
        if k ∉ discared
            if k == "nest"
                d[:nest_coverage] = v
            else
                d[Symbol(k)] = v
            end
        end
    end
    return (; pairs(d)...)
end
data = deserialize("/home/yakir/coffeebeetlearticle/all/data")
data["nest#closed person#therese displace_direction#none displace_location#feeder"] = data["nest#closed person#therese"]
delete!(data, "nest#closed person#therese")
df = DataFrame(parsetitle(k, r) for (k, v) in data for r in v.runs)

switchdirections(_::Missing) = missing
switchdirections(d) =   d == "left" ? "right" :
                        d == "right" ? "left" :
                        # d == "towards" ? "away" :
                        # d == "away" ? "towards" : 
                        d
getgroup(displace_location::Missing, transfer, displace_direction) = transfer
getgroup(displace_location, transfer, displace_direction) = displace_location == "nest" ? "zero" : displace_direction
getset(_::Missing, d) = d == "none" ? "Closed" : "Displacement"
getset(_, __) = "Transfer"

@. df[!, :displace_direction] = switchdirections.(df.displace_direction)
@. df[!, :group] .= getgroup(df.displace_location, df.transfer, df.displace_direction)
@. df[!, :set] = getset(df.transfer, df.group)

categorical!(df, [:group, :set, :displace_direction, :displace_location, :nest_coverage, :transfer])
levels!(df.group, ["none", "left", "right", "away", "towards", "zero", "back", "far"])

filter!(r -> r.group ≠ "far" || r.title == "transfer#far person#therese", df)
df[!, :direction_deviation]  = [angle(r.fictive_nest - r.feeder, turningpoint(r.track) - r.feeder) for r in eachrow(df)]
max_direction_deviation = maximum(r.direction_deviation for r in eachrow(df) if r.group ∉ ("far", "zero"))
mean_direction_deviation = mean(r.direction_deviation for r in eachrow(df) if r.group ∉ ("far", "zero"))
filter!(r -> r.group ≠ "far" || r.direction_deviation < 4mean_direction_deviation, df)

intended(d::AbstractString) =  d == "none" ? DungAnalyse.Point(0,0) :
                d == "away" ? DungAnalyse.Point(0, -50) :
                d == "towards" ? DungAnalyse.Point(0, 50) :
                d == "right" ? DungAnalyse.Point(50, 0) :
                d == "left" ? DungAnalyse.Point(-50, 0) :
                d == "zero" ? DungAnalyse.Point(0, -130) :
                d == "back" ? DungAnalyse.Point(0,0) :
                d == "far" ? DungAnalyse.Point(0,0) :
                error("unknown displacement")
intended(d) = intended(string(d))

# df[!, :intended] .= intendedf.(df.group)
# categorical!(df, :intended)

_get_rotation_center(displace_location::Missing, nest, fictive_nest) = fictive_nest
_get_rotation_center(displace_location, nest, fictive_nest) = displace_location == "feeder" ? fictive_nest : nest
_get_zeroing(nest::Missing, fictive_nest) = fictive_nest
_get_zeroing(nest, fictive_nest) = nest
function createtrans(nest, displace_location, fictive_nest, feeder)
    v = feeder - _get_rotation_center(displace_location, nest, fictive_nest)
    α = atan(v[2], v[1])
    rot = LinearMap(Angle2d(-π/2 - α))
    trans = Translation(-_get_zeroing(nest, fictive_nest))
    passmissing(rot ∘ trans)
end

df[!, :turning_point] .= zero.(df.feeder)
df[!, :center_of_search] .= zero.(df.feeder)
for r in eachrow(df)
    trans = createtrans(r.nest, r.displace_location, r.fictive_nest, r.feeder)
    @. r.track.coords = trans(r.track.coords)
    @. r.track.rawcoords.xy .= trans(r.track.rawcoords.xy)
    r.feeder = trans(r.feeder)
    r.fictive_nest = trans(r.fictive_nest)
    r.nest = trans(r.nest)
    Δ = intended(r.group) - r.fictive_nest
    r.turning_point = turningpoint(r.track) + Δ
    r.center_of_search = searchcenter(r.track) + Δ
end

function highlight(c, i, n)
    h = HSL(c)
    HSL(h.h, h.s, i/(n + 1))
end
groups = levels(df.group)
nc = length(groups)
colors = OrderedDict(zip(groups, [colorant"black"; distinguishable_colors(nc - 1, [colorant"white", colorant"black"], dropseed = true)]))

gdf = groupby(df, [:group, :nest_coverage])

function getcolor(g)
    n = length(g)
    c = colors[g[1]]
    [highlight(c, i, n) for i in eachindex(g)]
end

DataFrames.transform!(gdf, :group => getcolor => :color)

############################## Descriptive stats ################################################

d = stack(df, [:turning_point, :center_of_search], variable_name = :point_type, value_name = :point_xy)
gd = groupby(d, [:point_type, :group, :nest_coverage])
meanstd(x) = (μ = mean(x); Ref(μ .± std(x, mean = μ)))
g = copy(combine(gd, :point_xy => meanstd => :point_xy, nrow))
g.id = repeat(1:nrow(g)÷2, outer = 2)
d = unstack(g, [:group, :nest_coverage, :nrow], :point_type, :point_xy)
sort!(d, [:nest_coverage, :group])
select!(d, [1,2,4,5,3])
tbl = d


myformat(_::Missing) = "-"
myformat(x::Measurement) = string(round(Int, x.val), "±", round(Int, x.err))
myformat(xs::AbstractVector{String}) = string("(", xs[1], ",", xs[2], ")")
myformat(xs::AbstractVector{Float64}) = string("(", round(xs[1], digits = 2), ",", round(xs[2], digits = 2), ")")
myformat(xs::AbstractVector{Measurement{T}}) where {T <: Real} = myformat(myformat.(xs))

open(joinpath("tables", "table1.txt"), "w") do io
    pretty_table(io, tbl, ["Group" "Nest coverage" "Turning point"                        "Gravity center"                           "n";
                           ""  ""           "μ ± σ" "μ ± σ" ""], 
                 hlines = [1,7],
                 alignment = [:l, :l, :c, :c, :r],
                 formatters = (v,i,j) -> 3 ≤ j ≤ 4  ? myformat(v) : v
                )
end

########### speed

speed(d) = (norm(p2.xy - p1.xy)/(p2.t - p1.t) for (p1, p2) in zip(d, lag(d, -1, default = d[end])) if p1.xy ≠ p2.xy)
function foo(d)
    length(d) < 10 && return missing
    s = speed(d)
    μ = mean(s)
    μ ± std(s, mean = μ)
end
df.homing_speed = [foo(track.rawcoords[1:track.tp]) for track in df.track]
df.search_speed = [foo(track.rawcoords[track.tp:end]) for track in df.track]
v = combine(groupby(df, :group), :homing_speed => mean ∘ skipmissing, :search_speed => mean ∘ skipmissing)
mean(skipmissing(vcat(df.homing_speed, df.search_speed)))

################### displaced closed nest stats

roundsig(x) = x ≤ 0.05 ? @sprintf("%.2g", x) : "ns"

gdf = groupby(df, [:group, :nest_coverage])
g = gdf[[(group = group, nest_coverage = "closed") for group in ("none", "right","left","towards", "away")]]

HypothesisTests.pvalue(t::StatsModels.TableRegressionModel) = GLM.coeftable(t).cols[4][2]

n = 10^6
tbls = map((:turning_point, :center_of_search)) do point
    data = combine(g, point => (x -> x .- Ref(mean(x))) => :centered)
    x1 = [r.centered for r in eachrow(data) if r.group == "none"]
    x2 = [r.centered for r in eachrow(data) if r.group ≠ "none"]
    tx = ApproximatePermutationTest(first.(x1), first.(x2), var, n)
    ty = ApproximatePermutationTest(last.(x1), last.(x2), var, n)
    σ = pvalue.([tx, ty], tail = :left)
    nσ = length(x1) + length(x2)
    @printf "Is the %s variance of the none group, (%i, %i) significantly smaller than the variancve of the displaced groups, (%i, %i)? P = (%.2g, %.2g) (n = %i)\n" replace(string(point), '_' => ' ') var(x1)... var(x2)... σ... nσ
    data = combine(g, point => (x -> first.(x)) => :x, point => (x -> last.(x)) => :y, :group => (x -> first.(intended.(x))) => :ix, :group => (x -> last.(intended.(x))) => :iy)
    tx = lm(@formula(ix ~ x), data)
    ty = lm(@formula(iy ~ y), data)
    nμ = nrow(data)
    @printf "Is the effect (%.2f, %.2f) of the displacement on the %s significant? P = (%.2g, %.2g) (n = %i)\n" coef(tx)[2] coef(ty)[2] replace(string(point), '_' => ' ') pvalue(tx) pvalue(ty) nμ
    combine(g, [:group, point] => ((g, x) -> myformat(roundsig.(pvalue.([
                                                                         OneSampleTTest(first.(x), first(intended(g[1]))), 
                                                                         OneSampleTTest(last.(x),   last(intended(g[1])))
                                                                        ])))) => point, nrow => :n)
end
tbl = hcat(tbls[1][!, Not(All(:n, :nest_coverage))], tbls[2][!, Not(All(:group, :nest_coverage))])
m = Matrix(tbl)
writedlm(joinpath("tables", "table2.csv"), m, ',')


############### transfer far stats

gdf = groupby(df, :group)
group = "far"
g = gdf[(group = group,)]

P = Dict(point => myformat(pvalue.([OneSampleTTest(first.(g[!, point]), first(intended(group))), OneSampleTTest(last.(g[!, point]),   last(intended(group)))])) for point in (:turning_point, :center_of_search))
nrow(g)


############### transfer back stats

gdf = groupby(df, :group)
group = "back"
g = gdf[(group = group,)]

P = Dict(point => pvalue(ApproximatePermutationTest(first.(g[!, point]), last.(g[!, point]), var, n), tail = :left) for point in (:turning_point, :center_of_search))


P = Dict(point => myformat(pvalue.([OneSampleTTest(first.(g[!, point]), first(intended(group))), OneSampleTTest(last.(g[!, point]),   last(intended(group)))])) for point in (:turning_point, :center_of_search))
nrow(g)


######################## common traits for the plots 
#
set_theme!(
    font = "noto sans", #"Arial", # 
    fontsize = 10,
    resolution = (493.228346, 500.0),
    linewidth = 1,
    strokewidth = 1px, 
    markersize = 3px, 
    rowgap = Fixed(0), 
    titlegap = Fixed(5), 
    groupgap = Fixed(10), 
    titlehalign = :left, 
    gridshalign = :left,
    patchsize = (1, 1),
    LLegend = (markersize = 10px, markerstrokewidth = 1, patchsize = (10, 10), rowgap = Fixed(0), titlegap = Fixed(5), groupgap = Fixed(10), titlehalign = :left, gridshalign = :left, framecolor = :transparent), 
    LAxis = (xticklabelsize = 8, yticklabelsize = 8, xlabel = "X (cm)", ylabel = "Y (cm)", aspect = DataAspect())
)



# take out the open away

######################## closed nest plots #######

function unused_space(ax)
    used_space = ax.scene.px_area[]
    suggestedbbox = ax.layoutobservables.suggestedbbox[]
    difference = widths(suggestedbbox) .- widths(used_space)
    # Vec2f0(0, difference[2])
end
mydecompose(origin, radii) = [origin + radii .* Iterators.reverse(sincos(t)) for t in range(0, stop = 2π, length = 51)]
mydecompose(x) = mydecompose(x.origin, x.radii)
brighten(c, p = 0.5) = weighted_color_mean(p, c, colorant"white")
darken(c, p = 0.5) = weighted_color_mean(p, c, colorant"black")


legendmarkers = OrderedDict(
                            "burrow" => (color = :black, marker = '⋆', strokecolor = :black, markerstrokewidth = 0.5, strokewidth = 0.5, markersize = 25px),
                            "fictive burrow" => (color = :white, marker = '⋆', strokecolor = :black, markerstrokewidth = 0.5, strokewidth = 0.5, markersize = 15px),
                            "turning point" => (color = :black, marker = markers["turning_point"], strokecolor = :transparent, markersize = 5px),
                            "center of search" => (color = :black, marker = markers["center_of_search"], strokecolor = :transparent, markersize = 5px),
                            "μ ± FWHM" => [(color = brighten(colorant"black", 0.75), strokecolor = :transparent, polypoints = mydecompose(Point2f0(0.5, 0.5), Vec2f0(0.75, 0.5))),
                                           (color = :white, marker = '+', strokecolor = :transparent, markersize = 10px), 
                                          ])

gdf = groupby(df, :group)
g = DataFrame(gdf[(group = "none",)], copycols = false)
sort!(g, :turning_point, by = norm)
g.color .= getcolor(g.group)

function binit(track, h, nbins, m, M)
    o = Union{Variance, Missing}[Variance() for _ in 1:nbins]
    d = track.rawcoords
    to = findfirst(x -> x.xy[2] > 0, d) - 1 # -3 is good
    for (p1, p2) in Iterators.take(zip(d, lag(d, -1, default = d[to])), to)
        y = -(p2.xy[2] + p1.xy[2])/2
        if m < y < M
            i = StatsBase.binindex(h, y)
            v = norm(p2.xy - p1.xy)/(p2.t - p1.t)
            fit!(o[i], v)
        end
    end
    replace!(x -> nobs(x) < 2 ? missing : x, o)
end

nbins = 6
m, M = (0, 130)
bins = range(m, stop = M, length = nbins + 1)
h = Histogram(bins)
g.id = 1:nrow(g)
DataFrames.transform!(g, :id, :track => ByRow(x -> binit(x, h, nbins, m, M)) => :yv)

μ = [Variance() for _ in 1:nbins]
for i in 1:nbins
    reduce(merge!, skipmissing(yv[i] for yv in g.yv), init = μ[i])
end


bandcolor = RGB(only(distinguishable_colors(1, [colorant"white"; g.color], dropseed = true))) #:yellow
scene, layout = layoutscene(0)
ax = layout[1,1] = LAxis(scene, 
                         aspect = nothing, 
                         xlabel = "Distance to nest (cm)",
                         ylabel = "Speed (cm/s)",
                        )
mbins = StatsBase.midpoints(bins)
bh = band!(ax, mbins, mean.(μ) .- std.(μ), mean.(μ) .+ std.(μ), color = RGBA(bandcolor, 0.5))
lh = lines!(ax, mbins, mean.(μ), color = :white, linewidth = 2)
for r in eachrow(g)
    xy = [Point2f0(x, mean(y)) for (x,y) in zip(mbins, r.yv) if !ismissing(y)]
    scatterlines!(ax, xy, color = r.color, markersize = 6px, marker = '●')
end
ylims!(ax, 0, ax.limits[].origin[2] + ax.limits[].widths[2])
layout[1,2] = LLegend(scene, [[MarkerElement(color = colors["none"], marker = '●', strokecolor = :black, markerstrokewidth = 0, markersize = 6px), LineElement(linestyle = nothing, linewidth = 1, color = colors["none"])], [bh, lh]], ["individual", "μ±σ"], linewidth = 2, strokewidth = 1px, markersize = 3px, rowgap = Fixed(0), titlegap = Fixed(5), groupgap = Fixed(10), titlehalign = :left, gridshalign = :left);
FileIO.save("a.pdf", scene)
FileIO.save("figures/speed with $nbins bins.pdf", scene)


scene, layout = layoutscene(0)
ax = layout[1,1] = LAxis(scene)
for r in eachrow(g)
    lines!(ax, r.track.coords, color = r.color)
end
resize!(scene, round.(Int, size(scene) .- unused_space(ax))...)
FileIO.save("a.pdf", scene)
FileIO.save(joinpath("figures", "closed tracks.pdf"), scene)

d = g[1:3,:]
postTP = 25
scene, layout = layoutscene(0)
# bigone = GridLayout()
ax = layout[1:2,1] = LAxis(scene)
for r in eachrow(d)
    lines!(ax, r.track.coords, color = r.color)
end
rect = FRect2D(-10,-10,20,20)
spinecolor = RGB(only(distinguishable_colors(1, [colorant"white"; d.color], dropseed = true))) #:yellow
lines!(ax, rect, color = spinecolor)
hidexdecorations!(ax, ticklabels = false, ticks = false, grid = false)
ax = layout[2,2] = LAxis(scene, bottomspinecolor = spinecolor,
           topspinecolor = spinecolor,
           leftspinecolor = spinecolor,
           rightspinecolor = spinecolor)
for r in eachrow(d)
    lines!(ax, homing(r.track), color = r.color)
    lines!(ax, searching(r.track)[1:postTP], color = 1:postTP, colormap = [r.color, colorant"white"])
    scatter!(ax, [r.turning_point]; legendmarkers["turning point"]..., color = RGBA(r.color, 0.75))
end
scatter!(ax, [zero(Point2f0)]; legendmarkers["nest"]...)
ax.targetlimits[] = rect
hidexdecorations!(ax, ticklabels = false, ticks = false, grid = false)
hideydecorations!(ax, ticklabels = false, ticks = false, grid = false)
x = filter(kv -> first(kv) ∈ ("nest", "turning point"), legendmarkers)
x["track"] = (linestyle = nothing, color = :black)
layout[1,2] = LLegend(scene, apply_element(values(x)), collect(keys(x)))
layout[3,1:2] = LText(scene, "X (cm)")
# resize!(scene, round.(Int, size(scene) .- unused_space(ax))...)
FileIO.save("a.pdf", scene)
FileIO.save(joinpath("figures", "turning point zoom 3 closed tracks.pdf"), scene)


######################## plot ####################


apply_element(xs) = apply_element.(xs)
apply_element(x::NamedTuple) =  :marker ∈ keys(x) ? MarkerElement(; x...) :
                                :linestyle ∈ keys(x) ? LineElement(; x...) :
                                PolyElement(; x...)

gdf = groupby(df, [:group, :nest_coverage])
g = gdf[[(group = group, nest_coverage = "closed") for group in ("right","left","towards", "away")]]

polys = OrderedDict(string(k.group) => (color = colors[k.group], strokecolor = :transparent) for k in keys(g))


scene, layout = layoutscene(0)
ax = layout[1,1] = LAxis(scene)
for gg in g
    for r in eachrow(gg)
        lines!(ax, r.track.coords, color = r.color)
    end
end
leg = ([LineElement(linestyle = nothing, color = colors[k]) for k in getfield.(NamedTuple.(keys(g)), :group)], get.(getfield.(NamedTuple.(keys(g)), :group)))
layout[1, 2] = LLegend(scene, leg...)
layout[1, 1, TopLeft()] = LText(scene, "A");
FileIO.save(joinpath("figures", "closed displacement tracks.pdf"), scene)
FileIO.save("a.pdf", scene)

function getellipse(xy)
    n = length(xy)
    X = Array{Float64}(undef, 2, n)
    for i in 1:n
        X[:,i] = xy[i]
    end
    dis = fit(DiagNormal, X)
    radii = sqrt(2log(2))*sqrt.(var(dis)) # half the FWHM
    (origin = Point2f0(mean(dis)), radii = Vec2f0(radii))
end

scene, layout = layoutscene(0)
axs = []
for (i, point_type) in enumerate(("turning_point", "center_of_search"))
    ax = layout[1,i] = LAxis(scene)
    scatter!(ax, [zero(Point2f0)]; legendmarkers["nest"]...)
    for (k, gg) in pairs(g)
        xy = gg[!, point_type]
        ellipse = getellipse(xy)
        c = colors[k.group]
        poly!(ax, mydecompose(ellipse), color = brighten(c, 0.25))
        scatter!(ax, [ellipse.origin]; legendmarkers["μ ± FWHM"][2]...)
        scatter!(ax, xy; legendmarkers[replace(point_type, "_" => " ")]..., color = RGBA(c, 0.75))
        scatter!(ax, [Point2f0(intended(k.group))]; legendmarkers["fictive nest"]..., strokecolor = colors[k.group])
    end
    push!(axs, ax)
end
hideydecorations!(axs[2], grid = false)
hidexdecorations!.(axs, grid = false, ticklabels = false, ticks = false)
layout[2, 1:2] = LText(scene, "X (cm)");
linkaxes!(axs...)
label_c = layout[1, 1, TopLeft()] = LText(scene, "B");
label_d = layout[1, 2, TopLeft()] = LText(scene, "C");
layout[1, 3] = LLegend(scene, apply_element.(values.([polys, legendmarkers])), collect.(keys.([polys, legendmarkers])), ["Displacements", "Points"])
resize!(scene, round.(Int, size(scene) .- unused_space(axs[1]))...)
FileIO.save(joinpath("figures", "displacement turning point and gravity center.pdf"), scene)


# 3 open away tracks
g = DataFrame(gdf[(group = "away", nest_coverage = "open")], copycols = false)
function distance2nest(track)
    length(searching(track)) < 10 && return Inf
    t = homing(track)
    i = findfirst(>(0) ∘ last, t)
    isnothing(i) ? Inf : abs(first(t[i]))
end
sort!(g, :track, by = distance2nest)
d = g[1:3,:]
c = colors[d.group[1]]
scene, layout = layoutscene(0)
ax = layout[1:2,1] = LAxis(scene)
for r in eachrow(d)
    lines!(ax, r.track.coords, color = r.color)
    scatter!(ax, [turningpoint(r.track)]; legendmarkers["turning point"]..., color = RGBA(c, 0.75))
end
scatter!(ax, [zero(Point2f0)]; legendmarkers["nest"]...)
x = filter(kv -> first(kv) ∈ ("nest", "turning point"), legendmarkers)
x["track"] = (linestyle = nothing, color = :black)
layout[1, 2] = LLegend(scene, apply_element(values(x)), collect(keys(x)))
layout[1, 1, TopLeft()] = LText(scene, "D");
resize!(scene, round.(Int, size(scene) .- unused_space(ax))...)
FileIO.save(joinpath("figures", "3 open away tracks.pdf"), scene)

FileIO.save("a.pdf", scene)


##### zero vector


gdf = groupby(df, :group)
g = gdf[[("zero",)]]


scene, layout = layoutscene(0)
ax = layout[1,1] = LAxis(scene)
for gg in g
    for r in eachrow(gg)
        lines!(ax, r.track.coords, color = r.color)
    end
end
leg = ([LineElement(linestyle = nothing, color = colors[k]) for k in getfield.(NamedTuple.(keys(g)), :group)], get.(getfield.(NamedTuple.(keys(g)), :group)))
# layout[1, 2] = LLegend(scene, leg...)
layout[1, 1, TopLeft()] = LText(scene, "A");
FileIO.save(joinpath("figures", "zero tracks.pdf"), scene)
FileIO.save("a.pdf", scene)


scene, layout = layoutscene(0)
axs = []
for (i, point_type) in enumerate(("turning_point", "center_of_search"))
    ax = layout[1,i] = LAxis(scene)
    scatter!(ax, [zero(Point2f0)]; legendmarkers["nest"]...)
    for (k, gg) in pairs(g)
        xy = gg[!, point_type]
        ellipse = getellipse(xy)
        c = colors[k.group]
        poly!(ax, mydecompose(ellipse), color = brighten(c, 0.25))
        scatter!(ax, [ellipse.origin]; legendmarkers["μ ± FWHM"][2]...)
        scatter!(ax, xy; legendmarkers[replace(point_type, "_" => " ")]..., color = RGBA(c, 0.75))
        scatter!(ax, [Point2f0(intended(k.group))]; legendmarkers["fictive nest"]..., strokecolor = colors[k.group])
    end
    push!(axs, ax)
end
hideydecorations!(axs[2], grid = false)
hidexdecorations!.(axs, grid = false, ticklabels = false, ticks = false)
layout[2, 1:2] = LText(scene, "X (cm)");
linkaxes!(axs...)
label_c = layout[1, 1, TopLeft()] = LText(scene, "B");
label_d = layout[1, 2, TopLeft()] = LText(scene, "C");
layout[1, 3] = LLegend(scene, apply_element(values(legendmarkers)), collect(keys(legendmarkers)))
resize!(scene, round.(Int, size(scene) .- unused_space(axs[1]))...)
FileIO.save("a.pdf", scene)

FileIO.save(joinpath("figures", "zero turning point and gravity center.pdf"), scene)


##### far

g = gdf[[("far",)]]

scene, layout = layoutscene(0)
ax = layout[1,1] = LAxis(scene)
for gg in g
    for r in eachrow(gg)
        lines!(ax, r.track.coords, color = r.color)
    end
end
leg = ([LineElement(linestyle = nothing, color = colors[k]) for k in getfield.(NamedTuple.(keys(g)), :group)], get.(getfield.(NamedTuple.(keys(g)), :group)))
# layout[1, 2] = LLegend(scene, leg...)
layout[1, 1, TopLeft()] = LText(scene, "A");
FileIO.save(joinpath("figures", "far tracks.pdf"), scene)
FileIO.save("a.pdf", scene)


scene, layout = layoutscene(0)
axs = []
for (i, point_type) in enumerate(("turning_point", "center_of_search"))
    ax = layout[1,i] = LAxis(scene)
    scatter!(ax, [zero(Point2f0)]; legendmarkers["nest"]...)
    for (k, gg) in pairs(g)
        xy = gg[!, point_type]
        ellipse = getellipse(xy)
        c = colors[k.group]
        poly!(ax, mydecompose(ellipse), color = brighten(c, 0.25))
        scatter!(ax, [ellipse.origin]; legendmarkers["μ ± FWHM"][2]...)
        scatter!(ax, xy; legendmarkers[replace(point_type, "_" => " ")]..., color = RGBA(c, 0.75))
        scatter!(ax, [Point2f0(intended(k.group))]; legendmarkers["fictive nest"]..., strokecolor = colors[k.group])
    end
    push!(axs, ax)
end
hideydecorations!(axs[2], grid = false)
hidexdecorations!.(axs, grid = false, ticklabels = false, ticks = false)
layout[2, 1:2] = LText(scene, "X (cm)");
linkaxes!(axs...)
label_c = layout[1, 1, TopLeft()] = LText(scene, "B");
label_d = layout[1, 2, TopLeft()] = LText(scene, "C");
layout[1, 3] = LLegend(scene, apply_element(values(legendmarkers)), collect(keys(legendmarkers)))
resize!(scene, round.(Int, size(scene) .- unused_space(axs[1]))...)
FileIO.save("a.pdf", scene)
FileIO.save(joinpath("figures", "far turning point and gravity center.pdf"), scene)

























######################## Figure 5 ################
gdf = groupby(df, [:group, :nest_coverage])
g = gdf[[(group = group, nest_coverage = "closed") for group in ("right","left","towards", "away")]]
polys = OrderedDict(string(k.group) => (color = colors[k.group], strokecolor = :transparent) for k in keys(g))
scene, layout = layoutscene(0)
ax = layout[1,1] = LAxis(scene)
for gg in g
    for r in eachrow(gg)
        lines!(ax, r.track.coords, color = r.color, linewidth = 0.3)
    end
end
hidexdecorations!(ax, grid = false, ticklabels = false, ticks = false)
hideydecorations!(ax, grid = false, ticklabels = false, ticks = false)
# leg = ([LineElement(linestyle = nothing, color = colors[k]) for k in getfield.(NamedTuple.(keys(g)), :group)], get.(getfield.(NamedTuple.(keys(g)), :group)))
# layout[1, 2] = LLegend(scene, leg...)
layout[1, 1, TopLeft()] = LText(scene, "A");
axs = []
for (i, point_type) in enumerate(("turning_point", "center_of_search"))
    ax = layout[2,i] = LAxis(scene)
    scatter!(ax, [zero(Point2f0)]; legendmarkers["burrow"]...)
    for (k, gg) in pairs(g)
        xy = gg[!, point_type]
        ellipse = getellipse(xy)
        c = colors[k.group]
        poly!(ax, mydecompose(ellipse), color = brighten(c, 0.25))
        scatter!(ax, [ellipse.origin]; legendmarkers["μ ± FWHM"][2]...)
        scatter!(ax, xy; legendmarkers[replace(point_type, "_" => " ")]..., color = RGBA(c, 0.75))
        scatter!(ax, [Point2f0(intended(k.group))]; legendmarkers["fictive burrow"]..., strokecolor = colors[k.group])
    end
    push!(axs, ax)
end
hideydecorations!(axs[1], grid = false, ticklabels = false, ticks = false)
hideydecorations!(axs[2], grid = false)
hidexdecorations!.(axs, grid = false, ticklabels = false, ticks = false)
layout[3, 1:2] = LText(scene, "X (cm)");
linkaxes!(axs...)
label_c = layout[2, 1, TopLeft()] = LText(scene, "C");
label_d = layout[2, 2, TopLeft()] = LText(scene, "D");
g = DataFrame(gdf[(group = "away", nest_coverage = "open")], copycols = false)
sort!(g, :track, by = distance2nest)
d = g[1:3,:]
c = colors[d.group[1]]
dd = GridLayout(scene)
ax = dd[1,1] = LAxis(scene)
for r in eachrow(d)
    lines!(ax, r.track.coords, color = r.color, linewidth = 0.3)
    scatter!(ax, [turningpoint(r.track)]; legendmarkers["turning point"]..., color = RGBA(c, 0.75))
end
scatter!(ax, [zero(Point2f0)]; legendmarkers["burrow"]...)
hidexdecorations!(ax, grid = false, ticklabels = false, ticks = false)
hideydecorations!(ax, grid = false, ticklabels = false, ticks = false)
x = copy(legendmarkers);
x["track"] = (linestyle = nothing, color = :black);
dd[1, 2] = LLegend(scene, apply_element.(values.([polys, legendmarkers])), collect.(keys.([polys, legendmarkers])), ["Displacements", "Points"]);
layout[1,2] = dd;
layout[1, 2, TopLeft()] = LText(scene, "B");
layout[1:2, 0] = LText(scene, "Y (cm)", rotation = π/2);
# resize!(scene, round.(Int, size(scene) .- unused_space(ax))...)
FileIO.save("a.pdf", scene)
FileIO.save(joinpath("figures", "figure 5.pdf"), scene)





######################## Figure 6 ################
gdf = groupby(df, [:group])
g = gdf[[("zero",)]]
polys = OrderedDict(string(k.group) => (color = colors[k.group], strokecolor = :transparent) for k in keys(g))
scene, layout = layoutscene(0)
ax = layout[1,1] = LAxis(scene)
for gg in g
    for r in eachrow(gg)
        lines!(ax, r.track.coords, color = r.color, linewidth = 0.3)
    end
end
hidexdecorations!(ax, grid = false, ticklabels = false, ticks = false)
layout[1, 1, TopLeft()] = LText(scene, "A");
axs = []
for (i, point_type) in enumerate(("turning_point", "center_of_search"))
    ax = layout[1,i + 1] = LAxis(scene)
    scatter!(ax, [zero(Point2f0)]; legendmarkers["burrow"]...)
    for (k, gg) in pairs(g)
        xy = gg[!, point_type]
        ellipse = getellipse(xy)
        c = colors[k.group]
        poly!(ax, mydecompose(ellipse), color = brighten(c, 0.25))
        scatter!(ax, [ellipse.origin]; legendmarkers["μ ± FWHM"][2]...)
        scatter!(ax, xy; legendmarkers[replace(point_type, "_" => " ")]..., color = RGBA(c, 0.75))
        scatter!(ax, [Point2f0(intended(k.group))]; legendmarkers["fictive burrow"]..., strokecolor = colors[k.group])
    end
    push!(axs, ax)
end
hideydecorations!(axs[1], grid = false, ticklabels = false, ticks = false)
hideydecorations!(axs[2], grid = false)
hidexdecorations!.(axs, grid = false, ticklabels = false, ticks = false)
layout[2, 1:3] = LText(scene, "X (cm)");
linkaxes!(axs...)
label_c = layout[1, 2, TopLeft()] = LText(scene, "B");
label_d = layout[1, 3, TopLeft()] = LText(scene, "C");
x = copy(legendmarkers);
x["track"] = (linestyle = nothing, color = :black);
layout[1,4] = LLegend(scene, apply_element(values(legendmarkers)), collect(keys(legendmarkers)));
resize!(scene, round.(Int, size(scene) .- unused_space(ax))...)
FileIO.save("a.pdf", scene)
FileIO.save(joinpath("figures", "figure 6.pdf"), scene)




######################## Figure 7 ################
gdf = groupby(df, [:group])
g = gdf[[("far",)]]
polys = OrderedDict(string(k.group) => (color = colors[k.group], strokecolor = :transparent) for k in keys(g))
scene, layout = layoutscene(0)
ax = layout[1,1] = LAxis(scene)
for gg in g
    for r in eachrow(gg)
        lines!(ax, r.track.coords, color = r.color, linewidth = 0.3)
    end
end
hidexdecorations!(ax, grid = false, ticklabels = false, ticks = false)
layout[1, 1, TopLeft()] = LText(scene, "A");
axs = []
for (i, point_type) in enumerate(("turning_point", "center_of_search"))
    ax = layout[1,i + 1] = LAxis(scene)
    scatter!(ax, [zero(Point2f0)]; legendmarkers["burrow"]...)
    for (k, gg) in pairs(g)
        xy = gg[!, point_type]
        ellipse = getellipse(xy)
        c = colors[k.group]
        poly!(ax, mydecompose(ellipse), color = brighten(c, 0.25))
        scatter!(ax, [ellipse.origin]; legendmarkers["μ ± FWHM"][2]...)
        scatter!(ax, xy; legendmarkers[replace(point_type, "_" => " ")]..., color = RGBA(c, 0.75))
        scatter!(ax, [Point2f0(intended(k.group))]; legendmarkers["fictive burrow"]..., strokecolor = colors[k.group])
    end
    push!(axs, ax)
end
hideydecorations!(axs[1], grid = false, ticklabels = false, ticks = false)
hideydecorations!(axs[2], grid = false)
hidexdecorations!.(axs, grid = false, ticklabels = false, ticks = false)
layout[2, 1:3] = LText(scene, "X (cm)");
linkaxes!(axs...)
label_c = layout[1, 2, TopLeft()] = LText(scene, "B");
label_d = layout[1, 3, TopLeft()] = LText(scene, "C");
x = copy(legendmarkers);
x["track"] = (linestyle = nothing, color = :black);
layout[1,4] = LLegend(scene, apply_element(values(legendmarkers)), collect(keys(legendmarkers)));
resize!(scene, round.(Int, size(scene) .- unused_space(ax))...)
FileIO.save("a.pdf", scene)
FileIO.save(joinpath("figures", "figure 7.pdf"), scene)








######################## Figure 5 ################ trying with out he open towards
gdf = groupby(df, [:group, :nest_coverage])
g = gdf[[(group = group, nest_coverage = "closed") for group in ("right","left","towards", "away")]]
polys = OrderedDict(string(k.group) => (color = colors[k.group], strokecolor = :transparent) for k in keys(g))
scene, layout = layoutscene(0)#, resolution = (493.228346, 200.0))
ax = layout[1,1] = LAxis(scene)
for gg in g
    for r in eachrow(gg)
        lines!(ax, r.track.coords, color = r.color, linewidth = 0.3)
    end
end
hidexdecorations!(ax, grid = false, ticklabels = false, ticks = false)
# hideydecorations!(ax, grid = false, ticklabels = false, ticks = false)
# leg = ([LineElement(linestyle = nothing, color = colors[k]) for k in getfield.(NamedTuple.(keys(g)), :group)], get.(getfield.(NamedTuple.(keys(g)), :group)))
# layout[1, 2] = LLegend(scene, leg...)
layout[1, 1, TopLeft()] = LText(scene, "A");
axs = []
for (i, point_type) in enumerate(("turning_point", "center_of_search"))
    ax = layout[1, i + 1] = LAxis(scene)
    scatter!(ax, [zero(Point2f0)]; legendmarkers["nest"]...)
    for (k, gg) in pairs(g)
        xy = gg[!, point_type]
        ellipse = getellipse(xy)
        c = colors[k.group]
        poly!(ax, mydecompose(ellipse), color = brighten(c, 0.25))
        scatter!(ax, [ellipse.origin]; legendmarkers["μ ± FWHM"][2]...)
        scatter!(ax, xy; legendmarkers[replace(point_type, "_" => " ")]..., color = RGBA(c, 0.75))
        scatter!(ax, [Point2f0(intended(k.group))]; legendmarkers["fictive nest"]..., strokecolor = colors[k.group])
    end
    push!(axs, ax)
end
hideydecorations!(axs[1], grid = false, ticklabels = false, ticks = false)
hideydecorations!(axs[2], grid = false)
hidexdecorations!.(axs, grid = false, ticklabels = false, ticks = false)
layout[2, 1:3] = LText(scene, "X (cm)");
linkaxes!(axs...)
label_c = layout[1, 2, TopLeft()] = LText(scene, "B");
label_d = layout[1, 3, TopLeft()] = LText(scene, "C");
x = copy(legendmarkers);
x["track"] = (linestyle = nothing, color = :black);
layout[3,1:3] = LLegend(scene, apply_element.(values.([polys, x])), collect.(keys.([polys, x])), ["Displacements", "Points"], orientation = :horizontal, nbanks = 2);
resize!(scene, round.(Int, size(scene) .- unused_space(axs[1]))...)
FileIO.save("a.pdf", scene)
FileIO.save("a.png", scene)
# FileIO.save(joinpath("figures", "figure 5.pdf"), scene)































































# TODO: 
# 1. the arrow for the feeder
# 2. legend


shapes = [MarkerElement(color = :black, marker = '⋆', strokecolor = :black, markerstrokewidth = 0.5, markersize = 35px),#, markerstrokewidth = 1, markersize = 20px), 
          MarkerElement(color = :white, marker = '⋆', strokecolor = :black, markerstrokewidth = 0.5, markersize = 15px),#, markerstrokewidth = 1, markersize = 20px), 
          MarkerElement(color = :black, marker = '●', strokecolor = :transparent, markersize = 5px),#, markerstrokewidth = 0),#, markerstrokewidth = 1, markersize = 10px), 
          [PolyElement(color = brighten(colorant"black", 0.75), strokecolor = :transparent, polypoints = mydecompose(Point2f0(0.5, 0.5), Vec2f0(0.75, 0.5))),
           MarkerElement(color = :white, marker = '+', strokecolor = :transparent, markersize = 10px), 
          ]]

for d in groupby(df, :set)
    for g in groupby(d, :nest_coverage)

        scene, layout = layoutscene(0, fontsize = 10, font = "noto sans", resolution = (493.228346, 400.0));
        ax = layout[1,1] = LAxis(scene, 
                                 xlabel = "X (cm)",
                                 ylabel = "Y (cm)",
                                 xticklabelsize = 8,
                                 yticklabelsize = 8,
                                 aspect = DataAspect())
        for r in eachrow(g)
            lines!(ax, r.track.coords, color = colors[r.group])
        end
        leg = ([LineElement(linestyle = nothing, color = colors[k]) for k in unique(g.group)], unique(g.group), string(g.set[1]))
        layout[1, 2] = LLegend(scene, leg..., markersize = 10px, markerstrokewidth = 1, patchsize = (10, 10), rowgap = Fixed(0), titlegap = Fixed(5), groupgap = Fixed(10), titlehalign = :left, gridshalign = :left)
        FileIO.save("a.pdf", scene)

        FileIO.save(joinpath("figures", string(g.nest_coverage[1], " ", g.set[1], " tracks.pdf")), scene)
    end
    dfstacked = stack(d, [:turning_point, :gravity_center], variable_name = :point_type, value_name = :point_xy)
    ellipses = by(dfstacked, [:point_type, :nest_coverage, :group, :set]) do g
        n = length(g.point_xy)
        X = Array{Float64}(undef, 2, n)
        for i in 1:n
            X[:,i] = g.point_xy[i]
        end
        dis = fit(DiagNormal, X)
        radii = sqrt(2log(2))*sqrt.(var(dis)) # half the FWHM
        ellipse = (origin = Point2f0(mean(dis)), radius = Vec2f0(radii))
        (ellipse = ellipse, )
    end
    for g in groupby(ellipses, [:point_type, :nest_coverage])
        scene, layout = layoutscene(0, fontsize = 10, font = "arial", resolution = (493.228346, 350.0));

        ax = layout[1,1] = LAxis(scene, 
                                 xlabel = "X (cm)",
                                 ylabel = "Y (cm)",
                                 xticklabelsize = 8,
                                 yticklabelsize = 8,
                                 aspect = DataAspect()
                                 # autolimitaspect = 1
                                )
        for g in groupby(g, :group)
            c = only(g.ellipse)
            xy = mydecompose(c.origin, c.radius)
            poly!(ax, xy, color = brighten(colors[g.group[1]], 0.25))
        end
        scatter!(ax, getfield.(g.ellipse, :origin), color = :white, marker = '+', markersize = 10px)
        scatter!(ax, [zero(Point2f0)], color = :black, strokecolor = :black, marker = '⋆', strokewidth = 0.5, markersize = 35px)
        scatter!(ax, [Point2f0(intended[k]) for k in g.group], color = :white, strokecolor = [colors[k] for k in g.group], marker = '⋆', strokewidth = 0.5, markersize = 15px)
        for g in groupby(g, :group)
            xy = [Point2f0(r[g.point_type[1]]) for r in eachrow(d) if r.group == g.group[1] && r.nest_coverage == g.nest_coverage[1]]
            scatter!(ax, xy, color = RGBA(colors[g.group[1]], 0.75), marker = '●', markersize = 5px)
        end
        polys = [PolyElement(color = colors[k], strokecolor = :transparent) for k in unique(g.group)]
        leg = ([polys, shapes], [unique(g.group), ["nest", "fictive nest", replace(string(g.point_type[1]), '_' => ' '), "μ ± FWHM"]], [string(g.set[1]), "Shapes"])
        layout[1, 2] = LLegend(scene, leg..., markersize = 10px, markerstrokewidth = 1, patchsize = (10, 10), rowgap = Fixed(0), titlegap = Fixed(5), groupgap = Fixed(10), titlehalign = :left, gridshalign = :left)

        #=textlayer = Scene(scene, ax.scene.px_area, camera = campixel!, raw = true)
        topright = lift(textlayer.px_area) do w
        xy = widths(w)
        xy .- 0.05max(xy...)
        end
        text!(textlayer, "feeder", position = topright, align = (:right, :top), textsize = 10, font = "noto sans")
        lines!(textlayer, @lift([$topright, $topright .- (30, 0)]))
        scatter!(textlayer, @lift([$topright]), marker = '►', markersize = 10px)=#

        # arrows!(textlayer, @lift(Point2f0[$topright .- (20,0)]), [Point2f0(20,0)], arrowsize = 10)
        # arrows!(textlayer, rand(Point2f0, 3), 100rand(Vec2f0, 3), arrowsize = 50)
        # arrows!(textlayer, @lift([$topright]), [Vec2f0(1,0)], color = :red)
        FileIO.save(joinpath("figures", string(g.nest_coverage[1], " ", g.point_type[1], " ", g.set[1], ".pdf")), scene)
        # FileIO.save(joinpath("figures", string(g.nest_coverage[1], " ", g.point_type[1], " ", g.set[1], ".png")), scene)
        # FileIO.save("a.pdf", scene)
    end
end


############################## Closed nest ################################################

############################## open away ################################################

d = filter(r -> r.group == "away" && r.nest_coverage == "open", df)
nosearch = <(5)∘length∘searching
fellin = count(nosearch, d.track)
n = nrow(d)
filter!(r -> !nosearch(r.track), d)

c = colors[d.group[1]]
scene, layout = layoutscene(0, fontsize = 10, font = "noto sans", resolution = (493.228346, 400.0));
ax = layout[1,1] = LAxis(scene, 
                         xlabel = "X (cm)",
                         ylabel = "Y (cm)",
                         xticklabelsize = 8,
                         yticklabelsize = 8,
                         aspect = DataAspect());
for r in eachrow(d)
    lines!(ax, homing(r.track), color = c)
    lines!(ax, searching(r.track)[1:postTP], color = 1:postTP, colormap = [c, colorant"white"])
    scatter!(ax, [turningpoint(r.track)], color = RGBA(c, 0.75), marker = '●', markersize = 5px)
end
scatter!(ax, [zero(Point2f0)], color = :black, strokecolor = :black, marker = '⋆', strokewidth = 0.5, markersize = 35px)
w = widths(ax.limits[])[1]
x = ax.limits[].origin[1]
xl = max(-x, x + w)
xlims!(ax, (-xl, xl))
shapes = [MarkerElement(color = :black, marker = '⋆', strokecolor = :black, markerstrokewidth = 0.5, markersize = 35px),
          LineElement(linestyle = nothing, color = c),
          MarkerElement(color = RGBA(c, 0.75), marker = '●', strokecolor = :transparent, markersize = 5px)
         ];
leg = (shapes, ["nest", "track", "turning point"]);
layout[1, 2] = LLegend(scene, leg..., markersize = 10px, markerstrokewidth = 1, patchsize = (10, 10), rowgap = Fixed(0), titlegap = Fixed(5), groupgap = Fixed(10), titlehalign = :left, gridshalign = :left);
FileIO.save(joinpath("figures", "open away tracks.pdf"), scene)
FileIO.save("a.pdf", scene)



############################## transfer back ################################################
#TODO:

d = filter(r -> r.group == "back", df)
Δ = d.nest .- d.fictive_nest
"The fictive nest was $(myformat(-μ .± σ)) away from the nest"

d = filter(r -> r.group ∈ ("back", "far"), df)


n = 10^6
tbls = map((:turning_point, :gravity_center)) do point
    data = by(d, :group) do g
        # μ = mean(g[!, point])
        (centered = [p for p in g[!, point]], )
    end
    x1 = [r.centered for r in eachrow(data) if r.group == "back"]
    x2 = [r.centered for r in eachrow(data) if r.group == "far"]
    tx = ApproximatePermutationTest(first.(x1), first.(x2), var, n)
    ty = ApproximatePermutationTest(last.(x1), last.(x2), var, n)
    σ = pvalue.([tx, ty], tail = :left)
    nσ = length(x1) + length(x2)
    @printf "Is the %s variance of the back group, (%i, %i) significantly smaller than the variancve of the far group, (%i, %i)? P = (%.2g, %.2g) (n = %i)\n" replace(string(point), '_' => ' ') var(x1)... var(x2)... σ... nσ
    by(d, :group) do g
        tx = OneSampleTTest(first.(g[!, point]), first(intended[g.group[1]]))
        ty = OneSampleTTest(last.(g[!, point]), last(intended[g.group[1]]))
        (turning_point = myformat([roundsig(pvalue(tx)), roundsig(pvalue(ty))]), n = nrow(g))
    end
end

tbl = insertcols!(tbls[1], 3, :gravity_center => tbls[2][!,:turning_point])

m = vcat(reshape(names(tbls[1]), 1, :), Matrix(tbl))
writedlm(joinpath("tables", "table3.csv"), m, ',')




d = filter(r -> r.group ∈ ("none", "far"), df)
max_direction_deviation = maximum(r.direction_deviation for r in eachrow(df) if r.group ∉ ("far", "zero"))
mean_direction_deviation = mean(r.direction_deviation for r in eachrow(df) if r.group ∉ ("far", "zero"))
filter!(r -> r.direction_deviation < 4mean_direction_deviation, d)


n = 10^6
tbls = map((:turning_point, :gravity_center)) do point
    data = by(d, :group) do g
        # μ = mean(g[!, point])
        (centered = [p for p in g[!, point]], )
    end
    x1 = [r.centered for r in eachrow(data) if r.group == "none"]
    x2 = [r.centered for r in eachrow(data) if r.group == "far"]
    tx = ApproximatePermutationTest(first.(x1), first.(x2), var, n)
    ty = ApproximatePermutationTest(last.(x1), last.(x2), var, n)
    σ = pvalue.([tx, ty], tail = :left)
    nσ = length(x1) + length(x2)
    @printf "Is the %s variance of the none group, (%i, %i) significantly smaller than the variancve of the far group, (%i, %i)? P = (%.2g, %.2g) (n = %i)\n" replace(string(point), '_' => ' ') var(x1)... var(x2)... σ... nσ
    by(d, :group) do g
        tx = OneSampleTTest(first.(g[!, point]), first(intended[g.group[1]]))
        ty = OneSampleTTest(last.(g[!, point]), last(intended[g.group[1]]))
        (turning_point = myformat([roundsig(pvalue(tx)), roundsig(pvalue(ty))]), n = nrow(g))
    end
end




##################### check for handedness

# handedness:

function turnleft(t)
    n = 15
    tp = turningpoint(t)
    s = searching(t)
    if length(s) < n
        missing
    else
        sp = s[n]
        last(tp) > last(sp)
    end
end
# d = filter(r -> r.group ∈ ("left", "right", "away", "towards", "none"), df)
d = copy(df)
d[!, :turnleft] .= turnleft.(d.track)
by(d, :group) do g
    x = skipmissing(g.turnleft)
    count(x)/length(collect(x))
end

# turn to nest

nestonleft(_::Missing, __) = missing
nestonleft(nest, track) = last(nest) < last(turningpoint(track))
d[!, :nestonleft] .= nestonleft.(d.fictive_nest, d.track)
# d[!, :nestonleft] .= nestonleft.(d.nest, d.track)
turned2nest(_::Missing, __) = missing
turned2nest(_, __::Missing) = missing
function turned2nest(nestonleft, turnleft)
    if nestonleft
        turnleft
    else
        !turnleft
    end
end
by(d, :group) do g
    x = skipmissing(turned2nest.(g.nestonleft, g.turnleft))
    count(x)/length(collect(x))
end
x = skipmissing(turned2nest.(d.nestonleft, d.turnleft))
count(x)/length(collect(x))










##################### plot tracks with turning points

d = filter(r -> r.displace_direction == "towards", dropmissing(df, :displace_direction))
y = map(eachrow(d)) do r
    i = findfirst(p -> p[1] < 0, r.track.coords)
    if !isnothing(i)
        abs(r.track.coords[i][2])
    else
        missing
    end
end
i = sortperm(y)
dd = d[i,:]

#=scene = Scene(scale_plot = false)
for i in 1:2
# lines!(scene, homing(dd.track[i]))
lines!(scene, dd.track[i].coords[1:dd.track[i].tp + 25])
scatter!(scene, [dd.track[i].coords[dd.track[i].tp]], markersize = 1, color = :red)
end
ylims!(scene, -10, 10)
FileIO.save("a.pdf", scene)=#



d = filter(r -> r.displace_direction == "left" && r.nest_coverage == "open", dropmissing(df, :displace_direction))
d.turning_point

#=

using DataFrames, Distributions
ngrp = 2
n = 10
grps = string.(range('a', length = ngrp))
df = DataFrame(grp = String[], x = Float64[])
d = Normal(0, 0.1)
for _ in 1:n
push!(df, (grp = "zero", x = rand(d)))
end
for (i, grp) in enumerate(grps)
d = Normal(10^i, exp(i))
for _ in 1:n
push!(df, (grp = grp, x = rand(d)))
end
end

x1 = filter(r -> r.grp == "zero", df).x
x2 = filter(r -> r.grp == "zero", df).x

t = ApproximatePermutationTest(x1, x2, var, 10^5)
pvalue(t)=#
