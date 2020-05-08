# add Serialization, DataStructures, CoordinateTransformations, Rotations, DataFrames, Missings, Distributions, AngleBetweenVectors, LinearAlgebra, FileIO, Colors
# add https://github.com/yakir12/DungAnalyse.jl
# add AbstractPlotting#use-vertical-dims-from-font CairoMakie#jk/scale_text MakieLayout#master

# 8 -> 1:58
# 4 -> 1:52
# 2 -> 2:49

using DungAnalyse, Serialization, DataStructures, CoordinateTransformations, Rotations, DataFrames, Missings, Distributions, AngleBetweenVectors, LinearAlgebra
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
function switchdirections(d)
    if d == "left"
        "right"
    elseif d == "right"
        "left"
    elseif d == "towards"
        "away"
    elseif d == "away"
        "towards"
    else
        d
    end
end
df[!, :displace_direction] .= switchdirections.(df.displace_direction)

getgroup(r) = ismissing(r.displace_location) ? r.transfer : r.displace_location == "nest" ? "zero" : r.displace_direction
df[!, :group] .= getgroup.(eachrow(df))

_set(_::Missing, d) = d == "none" ? "Closed" : "Displacement"
_set(_, __) = "Transfer"
@. df[!, :set] = _set(df.transfer, df.group)

categorical!(df, [:group, :set, :displace_direction, :displace_location, :nest_coverage, :transfer])
levels!(df.group, ["none", "left", "right", "away", "towards", "zero", "back", "far"])
filter!(r -> r.group ≠ "far" || r.title == "transfer#far person#therese", df)

intended = Dict("none" => DungAnalyse.Point(0,0), 
                "away" => DungAnalyse.Point(0, 50),
                "towards" => DungAnalyse.Point(0, -50),
                "right" => DungAnalyse.Point(50, 0),
                "left" => DungAnalyse.Point(-50, 0),
                "zero" => DungAnalyse.Point(0, -130),
                "back" => DungAnalyse.Point(0,0),
                "far" => DungAnalyse.Point(0,0))

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

for r in eachrow(df)
    trans = createtrans(r.nest, r.displace_location, r.fictive_nest, r.feeder)
    @. r.track.coords = trans(r.track.coords)
    @. r.track.rawcoords.xy .= trans(r.track.rawcoords.xy)
    r.feeder = trans(r.feeder)
    r.fictive_nest = trans(r.fictive_nest)
    r.nest = trans(r.nest)
end

df[!, :turning_point] .= zero.(df.feeder)
df[!, :gravity_center] .= zero.(df.feeder)
df[!, :direction_deviation] .= 0.0
for r in eachrow(df)
    tp = turningpoint(r.track)
    r.turning_point = tp - r.fictive_nest + intended[r.group]
    r.gravity_center = searchcenter(r.track) - r.fictive_nest + intended[r.group]
    r.direction_deviation = angle(r.fictive_nest - r.feeder, tp - r.feeder)
end

max_direction_deviation = maximum(r.direction_deviation for r in eachrow(df) if r.group ∉ ("far", "zero"))
mean_direction_deviation = mean(r.direction_deviation for r in eachrow(df) if r.group ∉ ("far", "zero"))
# filter!(r -> r.direction_deviation < 4mean_direction_deviation, df)


############################## Descriptive stats ################################################

d = stack(df, [:turning_point, :gravity_center], variable_name = :point_type, value_name = :point_xy)
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

######################## common traits for the plots 
#
function highlight(c, i, n)
    h = HSL(c)
    HSL(h.h, h.s, i/(n + 1))
end
groups = levels(df.group)
nc = length(groups)
colors = OrderedDict(zip(groups, [colorant"black"; distinguishable_colors(nc - 1, [colorant"white", colorant"black"], dropseed = true)]))

gdf = groupby(df, [:group, :nest_coverage])
DataFrames.transform!(gdf, :group, nrow => :n)
DataFrames.transform!(gdf, :group => eachindex => :id)
DataFrames.transform!(gdf, :group => (g -> colors[g[1]]) => :groupcolor)
df.color = highlight.(df.groupcolor, df.id, df.n)

######################## closed nest plots #######

g = filter(r -> r.group == "none", df)
sort!(g, :turning_point, by = norm)

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
DataFrames.transform!(g, :id, :track => ByRow(x -> binit(x, h, nbins, m, M)) => :yv)

μ = [Variance() for _ in 1:nbins]
for i in 1:nbins
    reduce(merge!, skipmissing(yv[i] for yv in g.yv), init = μ[i])
end


bandcolor = RGB(only(distinguishable_colors(1, [colorant"white"; g.color], dropseed = true))) #:yellow
scene, layout = layoutscene(0, fontsize = 10, font = "noto sans", resolution = (493.228346, 400.0));
ax = layout[1,1] = LAxis(scene, 
                         xreversed = true,
                         xlabel = "Distance to nest (cm)",
                         ylabel = "Speed (cm/s)",
                         xticklabelsize = 8,
                         yticklabelsize = 8,
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


scene, layout = layoutscene(0, fontsize = 10, font = "noto sans", resolution = (493.228346, 400.0));
ax = layout[1,1] = LAxis(scene, 
                         xlabel = "X (cm)",
                         ylabel = "Y (cm)",
                         xticklabelsize = 8,
                         yticklabelsize = 8,
                         aspect = DataAspect())
for r in eachrow(g)
    lines!(ax, r.track.coords, color = r.color)
end
FileIO.save("a.pdf", scene)
FileIO.save(joinpath("figures", "closed tracks.pdf"), scene)

#=function distance2nest(track)
    t = homing(track)
    i = findfirst(>(0) ∘ last, t)
    abs(first(t[i]))
end=#
d = g[1:3,:]
postTP = 25
scene, layout = layoutscene(0, fontsize = 10, font = "noto sans", resolution = (493.228346, 493.228346));
# bigone = GridLayout()
ax = layout[1:2,1] = LAxis(scene, 
                         xlabel = "X (cm)",
                         ylabel = "Y (cm)",
                         xticklabelsize = 8,
                         yticklabelsize = 8,
                         aspect = DataAspect())
for r in eachrow(d)
    lines!(ax, r.track.coords, color = r.color)
end
rect = FRect2D(-10,-10,20,20)
spinecolor = RGB(only(distinguishable_colors(1, [colorant"white"; d.color], dropseed = true))) #:yellow
lines!(ax, rect, color = spinecolor)
ax = LAxis(scene, 
                         xlabel = "X (cm)",
                         xticklabelsize = 8,
                         yticklabelsize = 8,
                         bottomspinecolor = spinecolor,
                         topspinecolor = spinecolor,
                         leftspinecolor = spinecolor,
                         rightspinecolor = spinecolor,
                         aspect = DataAspect())
for r in eachrow(d)
    lines!(ax, homing(r.track), color = r.color)
    lines!(ax, searching(r.track)[1:postTP], color = 1:postTP, colormap = [r.color, colorant"white"])
    scatter!(ax, [r.turning_point], color = RGBA(r.color, 0.75), marker = '●', markersize = 5px)
end
scatter!(ax, [zero(Point2f0)], color = :black, strokecolor = :black, marker = '⋆', strokewidth = 0.5, markersize = 35px)
ax.targetlimits[] = rect
c = colors["none"]
shapes = [MarkerElement(color = :black, marker = '⋆', strokecolor = :black, markerstrokewidth = 0.5, markersize = 35px),
          LineElement(linestyle = nothing, color = c),
          MarkerElement(color = RGBA(c, 0.75), marker = '●', strokecolor = :transparent, markersize = 5px)
         ];
leg =  LLegend(scene, shapes, ["nest", "track", "turning point"], markersize = 10px, linewidth = 1, markerstrokewidth = 1, patchsize = (10, 10), rowgap = Fixed(0), titlegap = Fixed(5), groupgap = Fixed(10), titlehalign = :left, gridshalign = :left);
zoomed = GridLayout()
zoomed[1:2,1:2] = leg;
zoomed[2:3,1:2] = ax;
layout[1:2, 2] = zoomed;
FileIO.save("a.pdf", scene)
FileIO.save(joinpath("figures", "turning point zoom 3 closed tracks.pdf"), scene)



mydecompose(origin, radii) = [origin + radii .* Iterators.reverse(sincos(t)) for t in range(0, stop = 2π, length = 51)]
brighten(c, p = 0.5) = weighted_color_mean(p, c, colorant"white")
darken(c, p = 0.5) = weighted_color_mean(p, c, colorant"black")
markers = Dict("turning_point" => '●', "gravity_center" => '■')
shapes = [MarkerElement(color = :black, marker = '⋆', strokecolor = :black, markerstrokewidth = 0.5, markersize = 35px),#, markerstrokewidth = 1, markersize = 20px), 
          MarkerElement(color = :white, marker = '⋆', strokecolor = :black, markerstrokewidth = 0.5, markersize = 15px),#, markerstrokewidth = 1, markersize = 20px), 
          MarkerElement(color = :black, marker = markers["turning_point"], strokecolor = :transparent, markersize = 5px),#, markerstrokewidth = 0),#, markerstrokewidth = 1, markersize = 10px), 
          MarkerElement(color = :black, marker = markers["gravity_center"], strokecolor = :transparent, markersize = 5px),#, markerstrokewidth = 0),#, markerstrokewidth = 1, markersize = 10px), 
          [PolyElement(color = brighten(colorant"black", 0.75), strokecolor = :transparent, polypoints = mydecompose(Point2f0(0.5, 0.5), Vec2f0(0.75, 0.5))),
           MarkerElement(color = :white, marker = '+', strokecolor = :transparent, markersize = 10px), 
          ]]
d = g
dfstacked = stack(d, [:turning_point, :gravity_center], variable_name = :point_type, value_name = :point_xy)
ellipses = combine(groupby(dfstacked, :point_type)) do g
    n = length(g.point_xy)
    X = Array{Float64}(undef, 2, n)
    for i in 1:n
        X[:,i] = g.point_xy[i]
    end
    dis = fit(DiagNormal, X)
    radii = sqrt(2log(2))*sqrt.(var(dis)) # half the FWHM
    ellipse = (origin = Point2f0(mean(dis)), radius = Vec2f0(radii))
    (ellipse = ellipse, xy = Ref(g.point_xy))
end
scene, layout = layoutscene(0, fontsize = 10, font = "arial", resolution = (493.228346, 300.0));
axs = []
for (i, r) in enumerate(eachrow(ellipses))
    ax = layout[1,i] = LAxis(scene, 
                             # xlabel = "X (cm)",
                             # ylabel = "Y (cm)",
                             xticklabelsize = 8,
                             yticklabelsize = 8,
                             aspect = DataAspect()
                             # autolimitaspect = 1
                            )
    push!(axs, ax)
    c = r.ellipse
    xy = mydecompose(c.origin, c.radius)
    poly!(ax, xy, color = brighten(colors["none"], 0.25))
    scatter!(ax, [c.origin], color = :white, marker = '+', markersize = 10px)
    scatter!(ax, [zero(Point2f0)], color = :black, strokecolor = :black, marker = '⋆', strokewidth = 0.5, markersize = 35px)
    scatter!(ax, [Point2f0(intended["none"])], color = :white, strokecolor = colors["none"], marker = '⋆', strokewidth = 0.5, markersize = 15px)
    scatter!(ax, r.xy[], color = RGBA(colors["none"], 0.75), marker = markers[get(r.point_type)], markersize = 5px)
end
axs[1].ylabel = "Y (cm)"
hideydecorations!(axs[2], grid = false)
layout[2, 1:2] = LText(scene, "X (cm)");
linkaxes!(axs...)
label_c = layout[1, 1, TopLeft()] = LText(scene, "C");
label_d = layout[1, 2, TopLeft()] = LText(scene, "D");
polys = [PolyElement(color = colors[k], strokecolor = :transparent) for k in unique(g.group)];
leg = (shapes, ["nest", "fictive nest", "turning point", "gravity center", "μ ± FWHM"]);
layout[1, 3] = LLegend(scene, leg..., markersize = 10px, markerstrokewidth = 1, patchsize = (10, 10), rowgap = Fixed(0), titlegap = Fixed(5), groupgap = Fixed(10), titlehalign = :left, gridshalign = :left);
FileIO.save("a.pdf", scene)
FileIO.save(joinpath("figures", "none turning point and gravity center.pdf"), scene)

# possibel speed and direction plots
#=function midpoints(x)
n = length(x)
y = Vector{Float64}(undef, n - 1)
for i = 1:n-1
y[i] = (x[i] + x[i + 1])/2
end
return y
end

xtick = -40:20:40
xticklabels = string.(xtick)
xticklabels[3] = "TP"
meandegree(x) = rad2deg(angle(sum(exp, im*x)))
bins = -50:10:50
labels = string.(Int.(midpoints(bins)))

d = DataFrame(x = Int[], y = Float64[], id = Int[], color = RGB[])
categorical!(d, :id)
for (i, r) in enumerate(eachrow(g))
Δ = diff(r.track.coords)
l = norm.(Δ)
L = cumsum(l)
tp = r.track.tp
L .-= L[tp]
i1 = findfirst(L .> -50)
i2 = findlast(L .< 50)
L = L[i1:i2]
v = l[i1:i2]./step(r.track.t)
l = parse.(Int, collect(cut(L, bins, labels = labels)))
a = DataFrame(x = l, y = v)
mu = aggregate(a, :x, mean)
rename!(mu, :y_mean => :y)
# @pgf push!(p, Plot({color = c}, Table(mu)))
mu[!, :id] .= i
mu[!, :color] .= r.color
append!(d, mu)
end

d1 = by(d, :id) do g
(coords = [Point2f0.(zip(g.x, g.y))], color = [g.color[1]])
end

scene, layout = layoutscene(0, fontsize = 10, font = "noto sans", resolution = (493.228346, 400.0));
ax = layout[1,1] = LAxis(scene, 
xlabel = "X (cm)",
ylabel = "Y (cm)",
xticklabelsize = 8,
yticklabelsize = 8)
for r in eachrow(d1)
lines!(ax, r.coords, color = r.color)
end
FileIO.save("a.pdf", scene)
FileIO.save("a.png", scene, px_per_unit = 2)=#



######################## plot ####################

# d = filter(r -> r.displace_direction ≠ "none" && r.nest_coverage == "closed", dropmissing(df, :displace_direction))

d = groupby(df, [:displace_direction, :nest_coverage])
g = d[[(displace_direction = dd, nest_coverage = "closed") for dd in ("right","left","towards", "away")]]


scene, layout = layoutscene(0, fontsize = 10, font = "noto sans", resolution = (493.228346, 400.0));
ax = layout[1,1] = LAxis(scene, 
                         xlabel = "X (cm)",
                         ylabel = "Y (cm)",
                         xticklabelsize = 8,
                         yticklabelsize = 8,
                         aspect = DataAspect())

for gg in g
    for r in eachrow(gg)
        lines!(ax, r.track.coords, color = colors[r.group])
    end
end
leg = ([LineElement(linestyle = nothing, color = colors[k]) for k in getfield.(NamedTuple.(keys(g)), :displace_direction)], get.(getfield.(NamedTuple.(keys(g)), :displace_direction)))
layout[1, 2] = LLegend(scene, leg..., markersize = 10px, markerstrokewidth = 1, patchsize = (10, 10), rowgap = Fixed(0), titlegap = Fixed(5), groupgap = Fixed(10), titlehalign = :left, gridshalign = :left)
FileIO.save(joinpath("figures", "closed displacement tracks.pdf"), scene)




# FileIO.save("a.pdf", scene)



#=for r in eachrow(df)
trans = passmissing(Translation(intended[r.group] - r.fictive_nest))
@. r.track.coords = trans(r.track.coords)
r.feeder = trans(r.feeder)
r.fictive_nest = trans(r.fictive_nest)
r.nest = trans(r.nest)
end=#

# max_direction_deviation = maximum(r.direction_deviation for r in eachrow(df) if r.group ∉ ("far", "zero"))
# by(df, :group, :direction_deviation => x -> round(Int, rad2deg(maximum(x))))
# max_direction_deviation = maximum(r.direction_deviation for r in eachrow(df) if r.group == "back")

# x = filter(r -> r.group == "far", df)
# x[!, :direction_deviation] .= round.(Int, rad2deg.(x.direction_deviation))
# sort!(x, (:title, :direction_deviation))
# x[:, [:comment, :direction_deviation]]

# [r.metadata.setup[:azimuth] for r in data["transfer#far person#therese"].runs]

# scene, layout = layoutscene(0, fontsize = 10, font = "noto sans", resolution = (493.228346, 400.0));
# ax = layout[1, 1] = LAxis(scene, xlabel = "direction deviation (°)", ylabel = "TP to fictive nest (cm)")
# scatter!(ax, x.direction_deviation, norm.(turningpoint.(x.track)), markersize = 9px) 
# FileIO.save("a.pdf", scene)

# filter!(r -> r.group ≠ "far" || r.direction_deviation ≤ 2max_direction_deviation, df)

# d = dropmissing(df, :displace_direction)

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

roundsig(x) = x ≤ 0.05 ? @sprintf("%.2g", x) : "ns"

d = filter(r -> r.group ∉ ("far", "back", "zero") && r.nest_coverage == "closed", df)

HypothesisTests.pvalue(t::StatsModels.TableRegressionModel) = GLM.coeftable(t).cols[4][2]

n = 10^6
tbls = map((:turning_point, :gravity_center)) do point
    data = by(d, :group) do g
        μ = mean(g[!, point])
        (centered = [p - μ for p in g[!, point]], )
    end
    x1 = [r.centered for r in eachrow(data) if r.group == "none"]
    x2 = [r.centered for r in eachrow(data) if r.group ≠ "none"]
    tx = ApproximatePermutationTest(first.(x1), first.(x2), var, n)
    ty = ApproximatePermutationTest(last.(x1), last.(x2), var, n)
    σ = pvalue.([tx, ty], tail = :left)
    nσ = length(x1) + length(x2)
    @printf "Is the %s variance of the none group, (%i, %i) significantly smaller than the variancve of the displaced groups, (%i, %i)? P = (%.2g, %.2g) (n = %i)\n" replace(string(point), '_' => ' ') var(x1)... var(x2)... σ... nσ
    data = DataFrame((x = first(r[point]), y = first(intended[r.group])) for r in eachrow(d))
    tx = lm(@formula(y ~ x), data)
    data = DataFrame((x = last(r[point]), y = last(intended[r.group])) for r in eachrow(d))
    ty = lm(@formula(y ~ x), data)
    nμ = nrow(data)
    @printf "Is the effect (%.2f, %.2f) of the displacement on the %s significant? P = (%.2g, %.2g) (n = %i)\n" coef(tx)[2] coef(ty)[2] replace(string(point), '_' => ' ') pvalue(tx) pvalue(ty) nμ
    by(d, :group) do g
        tx = OneSampleTTest(first.(g[!, point]), first(intended[g.group[1]]))
        ty = OneSampleTTest(last.(g[!, point]), last(intended[g.group[1]]))
        (turning_point = myformat([roundsig(pvalue(tx)), roundsig(pvalue(ty))]), n = nrow(g))
    end
end
tbl = insertcols!(tbls[1], 3, :gravity_center => tbls[2][!,:turning_point])
m = vcat(reshape(names(tbls[1]), 1, :), Matrix(tbl))
writedlm(joinpath("tables", "table2.csv"), m, ',')

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
