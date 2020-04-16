using DungAnalyse, Serialization, DataStructures, CoordinateTransformations, Rotations, DataFrames, Missings, DataStructures, Distributions
using CairoMakie, MakieLayout, FileIO, DataStructures, AbstractPlotting, Colors
import AbstractPlotting:px
CairoMakie.activate!()
# include("plot_utils.jl")

keep = (:displace_location, :displace_direction , :nest_coverage, :transfer)
discared = ("person", "pellet")
function parsetitle(title, run)
    d = Dict{Symbol, Any}(k => missing for k in keep)
    d[:nest_coverage] = "open"
    d[:nest] = run.originalnest
    d[:feeder] = run.feeder
    d[:fictive_nest] = run.nest
    d[:track] = run.track
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
df = DataFrame(parsetitle(k, r.data) for (k, v) in data for r in v.runs)

_set(_::Missing) = "Displacement"
_set(_) = "Transfer"
@. df[!, :set] = _set(df.transfer)

getgroup(r) = ismissing(r.displace_location) ? r.transfer : r.displace_location == "nest" ? "zero" : r.displace_direction
df[!, :group] .= getgroup.(eachrow(df))

categorical!(df, [:group, :set, :displace_direction, :displace_location, :nest_coverage, :transfer])
levels!(df.group, ["none", "left", "right", "away", "towards", "zero", "back", "far"])

intended = Dict("none" => DungAnalyse.Point(0,0), 
                "towards" => DungAnalyse.Point(-50, 0),
                "away" => DungAnalyse.Point(50, 0),
                "left" => DungAnalyse.Point(0, 50),
                "right" => DungAnalyse.Point(0, -50),
                "zero" => DungAnalyse.Point(130, 0),
                "back" => DungAnalyse.Point(0,0),
                "far" => DungAnalyse.Point(0,0))


_get_rotation_center(displace_location, nest, fictive_nest) = fictive_nest
_get_rotation_center(displace_location::String, nest, fictive_nest) = displace_location == "feeder" ? fictive_nest : nest
_get_zeroing(nest::Missing, fictive_nest) = fictive_nest
_get_zeroing(nest, fictive_nest) = nest
function createtrans(nest, displace_location, fictive_nest, feeder)
    v = feeder - _get_rotation_center(displace_location, nest, fictive_nest)
    α = atan(v[2], v[1])
    rot = LinearMap(Angle2d(-α))
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

groups = levels(df.group)
nc = length(groups)
colors = OrderedDict(zip(groups, [colorant"gray"; distinguishable_colors(nc - 1, [colorant"white", colorant"black", colorant"gray"], dropseed = true)]))

for r in eachrow(df)
    trans = passmissing(Translation(intended[r.group] - r.fictive_nest))
    @. r.track.coords = trans(r.track.coords)
    r.feeder = trans(r.feeder)
    r.fictive_nest = trans(r.fictive_nest)
    r.nest = trans(r.nest)
end

df[!, :turning_point] .= turningpoint.(df.track)
df[!, :gravity_center] .= searchcenter.(df.track)

# d = dropmissing(df, :displace_direction)
mydecompose(origin, radii) = [origin + radii .* Iterators.reverse(sincos(t)) for t in range(0, stop = 2π, length = 51)]
brighten(c, p = 0.5) = weighted_color_mean(p, c, colorant"white")
darken(c, p = 0.5) = weighted_color_mean(p, c, colorant"black")

shapes = [MarkerElement(color = :black, marker = '▼', strokecolor = :black, markerstrokewidth = 1),#, markerstrokewidth = 1, markersize = 20px), 
          MarkerElement(color = brighten(colorant"black", 0.75), marker = '●', strokecolor = :black, markerstrokewidth = 1),#, markerstrokewidth = 1, markersize = 10px), 
          MarkerElement(color = brighten(colorant"black", 0.75), marker = '▲', strokecolor = :black, markerstrokewidth = 1),#, markerstrokewidth = 1, markersize = 10px), 
          [PolyElement(color = brighten(colorant"black", 0.75), strokecolor = :transparent, polypoints = mydecompose(Point2f0(0.5, 0.5), Vec2f0(0.5, 0.75))),
           MarkerElement(color = :white, marker = '+', strokecolor = :transparent, markersize = 10px), 
          ]]

for d in groupby(df, :set)
    for g in groupby(d, :nest_coverage)
        scene, layout = layoutscene(fontsize = 10, font = "noto sans", resolution = (493.228346, 400.0));
        ax = layout[1,1] = LAxis(scene, 
                                  xlabel = "X (cm)",
                                  ylabel = "Y (cm)",
                                  xticklabelsize = 8,
                                  yticklabelsize = 8,
                                  aspect = DataAspect())
        for r in eachrow(g)
            lines!(ax, r.track.coords, color = colors[r.group])
        end
        # xlims!(ax, -240, 240)
        # ylims!(ax, -200, 200)
        leg = ([LineElement(linestyle = nothing, color = colors[k]) for k in unique(g.group)], unique(g.group), string(g.set[1]))
        layout[1, 2] = LLegend(scene, leg..., markersize = 10px, markerstrokewidth = 1, patchsize = (10, 10), rowgap = Fixed(0), titlegap = Fixed(5), groupgap = Fixed(10), titlehalign = :left, gridshalign = :left)
        FileIO.save(joinpath("figures", string(g.nest_coverage[1], " ", g.set[1], " tracks.pdf")), scene)
        # FileIO.save("a.pdf", scene)
    end
    dfstacked = stack(d, [:turning_point, :gravity_center], variable_name = :point_type, value_name = :point_xy)
    ellipses = by(dfstacked, [:point_type, :nest_coverage, :group, :set]) do g
        n = length(g.point_xy)
        X = Array{Float64}(undef, 2, n)
        for i in 1:n
            X[:,i] = g.point_xy[i] - g.fictive_nest[i] + intended[g.group[1]]
        end
        dis = fit(DiagNormal, X)
        radii = sqrt(2log(2))*sqrt.(var(dis)) # half the FWHM
        ellipse = (origin = Point2f0(mean(dis)), radius = Vec2f0(radii))
        (ellipse = ellipse, )
    end
    for g in groupby(ellipses, [:point_type, :nest_coverage])
        scene, layout = layoutscene(fontsize = 10, font = "noto sans", resolution = (493.228346, 400.0));

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
            poly!(ax, xy, color = brighten(colors[g.group[1]]))
        end
        scatter!(ax, getfield.(g.ellipse, :origin), color = :white, marker = '+', markersize = 10px)
        scatter!(ax, [Point2f0(intended[k]) for k in g.group], color = [colors[k] for k in g.group], strokewidth = 1, strokecolor = [darken(colors[k]) for k in g.group], marker = '▲', markersize = 10px)
        for g in groupby(g, :group)
            xy = [Point2f0(r[g.point_type[1]]) for r in eachrow(d) if r.group == g.group[1] && r.nest_coverage == g.nest_coverage[1]]
            scatter!(ax, xy, color = colors[g.group[1]], marker = '●', strokewidth = 1, strokecolor = darken(colors[g.group[1]]), markersize = 10px)
        end
        scatter!(ax, [zero(Point2f0)], color = :black, marker = '▼', strokecolor = :black, strokewidth = 1, markersize = 10px)
        # xlims!(ax, -240, 240)
        # ylims!(ax, -200, 200)
        polys = [PolyElement(color = colors[k], strokecolor = :transparent) for k in unique(g.group)]
        leg = ([polys, shapes], [unique(g.group), ["nest", replace(string(g.point_type[1]), '_' => ' '), "fictive nest", "μ ± FWHM"]], [string(g.set[1]), "Shapes"])
        layout[1, 2] = LLegend(scene, leg..., markersize = 10px, markerstrokewidth = 1, patchsize = (10, 10), rowgap = Fixed(0), titlegap = Fixed(5), groupgap = Fixed(10), titlehalign = :left, gridshalign = :left)
        FileIO.save(joinpath("figures", string(g.nest_coverage[1], " ", g.point_type[1], " ", g.set[1], ".pdf")), scene)
        # FileIO.save(joinpath("figures", string(g.nest_coverage[1], " ", g.point_type[1], " ", g.set[1], ".png")), scene)
        # FileIO.save("a.pdf", scene)
    end
end



##################### plot tracks with turning points

