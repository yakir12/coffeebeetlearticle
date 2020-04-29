# add Serialization, DataStructures, CoordinateTransformations, Rotations, DataFrames, Missings, Distributions, AngleBetweenVectors, LinearAlgebra, FileIO, Colors
# add https://github.com/yakir12/DungAnalyse.jl
# add AbstractPlotting#use-vertical-dims-from-font CairoMakie#jk/scale_text MakieLayout#master

using DungAnalyse, Serialization, DataStructures, CoordinateTransformations, Rotations, DataFrames, Missings, Distributions, AngleBetweenVectors, LinearAlgebra
using CairoMakie, MakieLayout, FileIO, AbstractPlotting, Colors
import AbstractPlotting:px
CairoMakie.activate!()
# include("plot_utils.jl")

using PrettyTables, Measurements, HypothesisTests, GLM, MixedModels

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
_get_rotation_center(displace_location::AbstractString, nest, fictive_nest) = displace_location == "feeder" ? fictive_nest : nest
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

############################## Descriptive stats ################################################

d = stack(df, [:turning_point, :gravity_center], variable_name = :point_type, value_name = :point_xy)
sort!(d, [:point_type, :nest_coverage, :group])
tbls = []
for g in groupby(d, :point_type)
    tbl = by(g, [:group, :nest_coverage]) do g
        d = Dict(g.point_type[1] => [mean(g.point_xy) .± std(g.point_xy)], :n => [nrow(g)])
        (; d...)
    end
    push!(tbls, tbl)
end

tp = tbls[2]
gc = tbls[1]
tbl = hcat(tp[!, Not(All(:n))], gc[:, Not(All(:group, :nest_coverage))])

myformat(_::Missing) = "-"
myformat(x::Measurement) = string(round(Int, x.val), "±", round(Int, x.err))
myformat(xs::AbstractVector{String}) = string("(", xs[1], ",", xs[2], ")")
myformat(xs::AbstractVector{Measurement{T}}) where {T <: Real} = myformat(myformat.(xs))

open("table1.txt", "w") do io
    pretty_table(io, tbl, ["Group" "Nest coverage" "Turning point"                        "Gravity center"                           "n";
                           ""  ""           "μ ± σ" "μ ± σ" ""], 
                 hlines = [1,7],
                 alignment = [:l, :l, :c, :c, :r],
                 formatters = (v,i,j) -> 3 ≤ j ≤ 4  ? myformat(v) : v
                )
end

######################## common traits for the plots 
#
groups = levels(df.group)
nc = length(groups)
colors = OrderedDict(zip(groups, [colorant"gray"; distinguishable_colors(nc - 1, [colorant"white", colorant"black", colorant"gray"], dropseed = true)]))

######################## closed nest plots #######

g = filter(r -> r.group == "none", df)
sort!(g, :turning_point, by =  norm)
g[!, :color] .= map(1:nrow(g)) do i
    c = HSL(colors["none"])
    HSL(c.h, c.s, (i - 1)/nrow(g))
end

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
FileIO.save(joinpath("figures", "closed tracks.pdf"), scene)

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
shapes = [MarkerElement(color = :black, marker = '⋆', strokecolor = :black, markerstrokewidth = 0.5, markersize = 35px),
          LineElement(linestyle = nothing, color = c),
          MarkerElement(color = RGBA(c, 0.75), marker = '●', strokecolor = :transparent, markersize = 5px)
         ];
leg =  LLegend(scene, shapes, ["nest", "track", "turning point"], markersize = 10px, markerstrokewidth = 1, patchsize = (10, 10), rowgap = Fixed(0), titlegap = Fixed(5), groupgap = Fixed(10), titlehalign = :left, gridshalign = :left);
zoomed = GridLayout()
zoomed[1:2,1:2] = leg;
zoomed[2:3,1:2] = ax;
layout[1:2, 2] = zoomed;
FileIO.save(joinpath("figures", "turning point zoom 3 closed tracks.pdf"), scene)


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

#TODO: 
# 1. the arrow for the feeder
# 2. legend
#
mydecompose(origin, radii) = [origin + radii .* Iterators.reverse(sincos(t)) for t in range(0, stop = 2π, length = 51)]
brighten(c, p = 0.5) = weighted_color_mean(p, c, colorant"white")
darken(c, p = 0.5) = weighted_color_mean(p, c, colorant"black")

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

function roundsig(x)
    sigs = OrderedDict("***" => 0.001, "**" => 0.01, "*" => 0.05, "ns" => 1)
    for (k, v) in sigs
        x ≤ v && return k
    end
end

d = filter(r -> r.group ∉ ("far", "back", "zero") && r.nest_coverage == "closed", df)

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
    σ = pvalue.([tx, ty])
    nσ = length(x1) + length(x2)
    data = DataFrame((x = first(r.turning_point), y = first(intended[r.group])) for r in eachrow(d))
    tx = lm(@formula(y ~ x), data)
    data = DataFrame((x = last(r.turning_point), y = last(intended[r.group])) for r in eachrow(d))
    ty = lm(@formula(y ~ x), data)
    μ = [GLM.coeftable(tx).cols[4][2], GLM.coeftable(ty).cols[4][2]]
    nμ = nrow(data)
    fn = by(d, :group) do g
        tx = OneSampleTTest(first.(g[!, point]), first(intended[g.group[1]]))
        ty = OneSampleTTest(last.(g[!, point]), last(intended[g.group[1]]))
        (Px = pvalue(tx), Py = pvalue(ty), n = nrow(g))
    end
    tbl = DataFrame(axis = ["x", "y", "n"], σ = [roundsig.(σ); nσ], μ = [roundsig.(μ); nμ])
    for r in eachrow(fn)
        tbl[!, Symbol(r.group)] .= [roundsig.([r.Px, r.Py]); r.n]
    end
    tbl
end

tbl = vcat(reshape(["Turning point"; fill("", 7)], 1, :), 
     reshape(names(tbls[1]), 1, :),
     Matrix(tbls[1][1:end - 1, :]),
     reshape(["Gravity center"; fill("", 7)], 1, :), 
     Matrix(tbls[2])
    )

writedlm("table1.csv", tbl, ',')

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
# 1. figure out the error
# 2. do some stats?

d = filter(r -> r.group == "back", df)
for r in eachrow(d)
    tp = turningpoint(r.track)
    r.turning_point = tp - r.nest + intended[r.group]
    r.gravity_center = searchcenter(r.track) - r.nest + intended[r.group]
end
for point in (:turning_point, :gravity_center)
    tp = d[!, point]
    n = length(tp)
    μx, μy = mean(tp)
    σx, σy = std(tp)
    @printf "The %s was located %i ± %i cm %s the nest and %i ± %i cm to the %s of it (mean ± std; n = %i)\n" replace(string(point), '_' => ' ') abs(μy) σy (μy > 0 ? "after" : "before") abs(μx) σx (μx > 0 ? "right" : "left") n
end
Δ = d.fictive_nest .- d.nest
μ = mean(Δ)
σ = std(Δ)
"The fictive nest was $(myformat(μ .± σ)) away from the nest"

ϵ = μ .± σ

d = filter(r -> r.group == "far", df)
for point in (:turning_point, :gravity_center)
    tp = [x[point] for x in eachrow(d)]
    n = length(tp)
    μx, μy = mean(tp)
    σx, σy = std(tp)
    @printf "The %s was located %i ± %i cm %s the nest and %i ± %i cm to the %s of it (mean ± std; n = %i)\n" replace(string(point), '_' => ' ') abs(μy) σy (μy > 0 ? "after" : "before") abs(μx) σx (μx > 0 ? "right" : "left") n
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
