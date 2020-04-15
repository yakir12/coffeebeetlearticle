Stacktrace:
 [1] exec(::Ptr{Nothing}, ::String, ::Int64, ::UInt32, ::Ptr{Nothing}) at ./pcre.jl:155
 [2] exec_r_data at ./pcre.jl:172 [inlined]
 [3] match(::Regex, ::String, ::Int64, ::UInt32) at ./regex.jl:277
 [4] match at ./regex.jl:275 [inlined]
 [5] match at ./regex.jl:294 [inlined]
 [6] _parse_colorant(::String) at /home/yakir/.julia/packages/Colors/P2rUx/src/parse.jl:58
 [7] _parse_colorant(::Type{RGBA{Float32}}, ::Type{ColorAlpha{RGB{Float32},Float32,4}}, ::String) at /home/yakir/.julia/packages/Colors/P2rUx/src/parse.jl:144
 [8] parse at /home/yakir/.julia/packages/Colors/P2rUx/src/parse.jl:181 [inlined]
 [9] convert_attribute at /home/yakir/.julia/packages/AbstractPlotting/WwXMh/src/conversions.jl:562 [inlined]
 [10] convert_attribute at /home/yakir/.julia/packages/AbstractPlotting/WwXMh/src/conversions.jl:559 [inlined]
 [11] to_color at /home/yakir/.julia/packages/AbstractPlotting/WwXMh/src/conversions.jl:13 [inlined]
 [12] (::AbstractPlotting.var"#380#382"{StaticArrays.SArray{Tuple{4,4},Float32,2,16},Base.GenericIOBuffer{Array{UInt8,1}},Array{Point{3,Float32},1},Array{RGBA{Float32},1},Array{Vec{2,Float32},1},Array{FreeTypeAbstraction.FTFont,1},Array{Quaternion{Float32},1}})(::Int64, ::FreeTypeAbstraction.FTFont, ::Tuple{String,Point{2,Float32}}, ::Symbol, ::Int64, ::Tuple{Symbol,Symbol}, ::Float32) at /home/yakir/.julia/packages/AbstractPlotting/WwXMh/src/basic_recipes/basic_recipes.jl:425
 [13] broadcast_foreach(::Function, ::UnitRange{Int64}, ::Vararg{Any,N} where N) at /home/yakir/.julia/packages/AbstractPlotting/WwXMh/src/utilities/utilities.jl:156
 [14] (::AbstractPlotting.var"#379#381"{Array{Point{3,Float32},1},Array{RGBA{Float32},1},Array{Vec{2,Float32},1},Array{FreeTypeAbstraction.FTFont,1},Array{Quaternion{Float32},1},Text{...}})(::StaticArrays.SArray{Tuple{4,4},Float32,2,16}, ::String, ::Array{Tuple{String,Point{2,Float32}},1}, ::Symbol, ::Vararg{Any,N} where N) at /home/yakir/.julia/packages/AbstractPlotting/WwXMh/src/basic_recipes/basic_recipes.jl:424
 [15] (::Observables.OnUpdate{AbstractPlotting.var"#379#381"{Array{Point{3,Float32},1},Array{RGBA{Float32},1},Array{Vec{2,Float32},1},Array{FreeTypeAbstraction.FTFont,1},Array{Quaternion{Float32},1},Text{...}},Tuple{Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Array{Tuple{String,Point{2,Float32}},1}},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any}}})(::Array{Tuple{String,Point{2,Float32}},1}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [16] setindex!(::Observables.Observable{Array{Tuple{String,Point{2,Float32}},1}}, ::Array{Tuple{String,Point{2,Float32}},1}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [17] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [18] MapUpdater at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:241 [inlined]
 [19] (::Observables.OnUpdate{Observables.MapUpdater{AbstractPlotting.var"#187#189"{Int64},Array{Tuple{String,Point{2,Float32}},1}},Tuple{Observables.Observable{Tuple{Array{Tuple{String,Point{2,Float32}},1}}}}})(::Tuple{Array{Tuple{String,Point{2,Float32}},1}}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [20] setindex!(::Observables.Observable{Tuple{Array{Tuple{String,Point{2,Float32}},1}}}, ::Tuple{Array{Tuple{String,Point{2,Float32}},1}}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [21] setindex!(::Observables.Observable{Tuple{Array{Tuple{String,Point{2,Float32}},1}}}, ::Tuple{Array{Tuple{String,Point{2,Float32}},1}}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126
 [22] (::AbstractPlotting.var"#211#213"{DataType,Observables.Observable{Tuple{Array{Tuple{String,Point{2,Float32}},1}}}})(::Tuple{}, ::Tuple{Array{Tuple{String,Point{2,Float32}},1}}) at /home/yakir/.julia/packages/AbstractPlotting/WwXMh/src/interfaces.jl:561
 [23] (::Observables.OnUpdate{AbstractPlotting.var"#211#213"{DataType,Observables.Observable{Tuple{Array{Tuple{String,Point{2,Float32}},1}}}},Tuple{Observables.Observable{Tuple{}},Observables.Observable{Tuple{Array{Tuple{String,Point{2,Float32}},1}}}}})(::Tuple{Array{Tuple{String,Point{2,Float32}},1}}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [24] setindex!(::Observables.Observable{Tuple{Array{Tuple{String,Point{2,Float32}},1}}}, ::Tuple{Array{Tuple{String,Point{2,Float32}},1}}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [25] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [26] MapUpdater at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:241 [inlined]
 [27] (::Observables.OnUpdate{Observables.MapUpdater{typeof(tuple),Tuple{Array{Tuple{String,Point{2,Float32}},1}}},Tuple{Observables.Observable{Array{Tuple{String,Point{2,Float32}},1}}}})(::Array{Tuple{String,Point{2,Float32}},1}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [28] setindex!(::Observables.Observable{Array{Tuple{String,Point{2,Float32}},1}}, ::Array{Tuple{String,Point{2,Float32}},1}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [29] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [30] (::MakieLayout.var"#89#109"{Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Observables.Observable{Array{Tuple{String,Point{2,Float32}},1}},Observables.Observable{Float32},Observables.Observable{Array{Point{2,Float32},1}}})(::Array{String,1}, ::Float32, ::Bool) at /home/yakir/.julia/packages/MakieLayout/bYPho/src/lineaxis.jl:199
 [31] (::Observables.OnUpdate{MakieLayout.var"#89#109"{Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Observables.Observable{Array{Tuple{String,Point{2,Float32}},1}},Observables.Observable{Float32},Observables.Observable{Array{Point{2,Float32},1}}},Tuple{Observables.Observable{Array{String,1}},Observables.Observable{Float32},Observables.Observable{Any}}})(::Array{String,1}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [32] setindex!(::Observables.Observable{Array{String,1}}, ::Array{String,1}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [33] setindex!(::Observables.Observable{Array{String,1}}, ::Array{String,1}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126
 [34] (::MakieLayout.var"#86#106"{Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Observables.Observable{Array{Point{2,Float32},1}},Observables.Observable{Array{String,1}}})(::Array{Float64,1}) at /home/yakir/.julia/packages/MakieLayout/bYPho/src/lineaxis.jl:178
 [35] (::Observables.OnUpdate{MakieLayout.var"#86#106"{Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Observables.Observable{Array{Point{2,Float32},1}},Observables.Observable{Array{String,1}}},Tuple{Observables.Observable{Array{Float64,1}}}})(::Array{Float64,1}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [36] setindex!(::Observables.Observable{Array{Float64,1}}, ::Array{Float64,1}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [37] setindex!(::Observables.Observable{Array{Float64,1}}, ::Array{Float64,1}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126
 [38] (::Observables.MapUpdater{MakieLayout.var"#85#105",Array{Float64,1}})(::Tuple{Float32,Tuple{Float32,Float32},Bool}, ::Vararg{Any,N} where N) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:241
 [39] (::Observables.OnUpdate{Observables.MapUpdater{MakieLayout.var"#85#105",Array{Float64,1}},Tuple{Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Observables.Observable{Any},Observables.Observable{Any}}})(::Tuple{Float32,Tuple{Float32,Float32},Bool}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [40] setindex!(::Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}}, ::Tuple{Float32,Tuple{Float32,Float32},Bool}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [41] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 ... (the last 4 lines are repeated 1 more time)
 [46] (::Observables.var"#3#4"{Any,Observables.Observable{Any}})(::Tuple{Point{2,Float32},Point{2,Float32}}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:47
 [47] #invokelatest#1 at ./essentials.jl:712 [inlined]
 [48] invokelatest at ./essentials.jl:711 [inlined]
 [49] setindex!(::Observables.Observable{Tuple{Point{2,Float32},Point{2,Float32}}}, ::Tuple{Point{2,Float32},Point{2,Float32}}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:132
 [50] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [51] (::Observables.MapUpdater{MakieLayout.var"#121#144",Tuple{Point{2,Float32},Point{2,Float32}}})(::Symbol, ::Vararg{Any,N} where N) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:241
 [52] (::Observables.OnUpdate{Observables.MapUpdater{MakieLayout.var"#121#144",Tuple{Point{2,Float32},Point{2,Float32}}},Tuple{Observables.Observable{Any},Observables.Observable{GeometryBasics.HyperRectangle{2,Int64}}}})(::GeometryBasics.HyperRectangle{2,Int64}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [53] setindex!(::Observables.Observable{GeometryBasics.HyperRectangle{2,Int64}}, ::GeometryBasics.HyperRectangle{2,Int64}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [54] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [55] MapUpdater at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:241 [inlined]
 [56] (::Observables.OnUpdate{Observables.MapUpdater{AbstractPlotting.var"#67#68",GeometryBasics.HyperRectangle{2,Int64}},Tuple{Observables.Observable{GeometryBasics.HyperRectangle{2,Int64}},Observables.Observable{GeometryBasics.HyperRectangle{2,Int64}}}})(::GeometryBasics.HyperRectangle{2,Int64}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [57] setindex!(::Observables.Observable{GeometryBasics.HyperRectangle{2,Int64}}, ::GeometryBasics.HyperRectangle{2,Int64}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [58] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [59] (::MakieLayout.var"#1#2"{Observables.Observable{GeometryBasics.HyperRectangle{2,Int64}}})(::GeometryBasics.HyperRectangle{2,Float32}, ::GeometryBasics.HyperRectangle{2,Float32}, ::Nothing) at /home/yakir/.julia/packages/MakieLayout/bYPho/src/helpers.jl:52
 [60] (::Observables.OnUpdate{MakieLayout.var"#1#2"{Observables.Observable{GeometryBasics.HyperRectangle{2,Int64}}},Tuple{Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}},Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}},Observables.Observable{Any}}})(::GeometryBasics.HyperRectangle{2,Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [61] setindex!(::Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}, ::GeometryBasics.HyperRectangle{2,Float32}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [62] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [63] (::GridLayoutBase.var"#115#116"{Observables.Observable{Tuple{Any,Any}},Observables.Observable{Tuple{Union{Nothing, Float32},Union{Nothing, Float32}}},Observables.Observable{Any},Observables.Observable{GridLayoutBase.RectSides{Float32}},Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}})(::GeometryBasics.HyperRectangle{2,Float32}, ::Tuple{Symbol,Symbol}, ::Tuple{Nothing,Nothing}) at /home/yakir/.julia/packages/GridLayoutBase/Y582M/src/layoutobservables.jl:267
 [64] (::Observables.OnUpdate{GridLayoutBase.var"#115#116"{Observables.Observable{Tuple{Any,Any}},Observables.Observable{Tuple{Union{Nothing, Float32},Union{Nothing, Float32}}},Observables.Observable{Any},Observables.Observable{GridLayoutBase.RectSides{Float32}},Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}},Tuple{Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}},Observables.Observable{Any},Observables.Observable{Tuple{Union{Nothing, Float32},Union{Nothing, Float32}}}}})(::GeometryBasics.HyperRectangle{2,Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [65] setindex!(::Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}, ::GeometryBasics.HyperRectangle{2,Float32}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [66] setindex!(::Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}, ::GeometryBasics.HyperRectangle{2,Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126
 [67] align_to_bbox!(::GridLayout, ::GeometryBasics.HyperRectangle{2,Float32}) at /home/yakir/.julia/packages/GridLayoutBase/Y582M/src/gridlayout.jl:682
 [68] (::GridLayoutBase.var"#6#8"{GridLayout})(::GeometryBasics.HyperRectangle{2,Float32}) at /home/yakir/.julia/packages/GridLayoutBase/Y582M/src/gridlayout.jl:36
 [69] #invokelatest#1 at ./essentials.jl:712 [inlined]
 [70] invokelatest at ./essentials.jl:711 [inlined]
 [71] setindex!(::Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}, ::GeometryBasics.HyperRectangle{2,Float32}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:132
 [72] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [73] (::GridLayoutBase.var"#115#116"{Observables.Observable{Tuple{Any,Any}},Observables.Observable{Tuple{Union{Nothing, Float32},Union{Nothing, Float32}}},Observables.Observable{GridLayoutBase.AlignMode},Observables.Observable{GridLayoutBase.RectSides{Float32}},Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}})(::GeometryBasics.HyperRectangle{2,Float32}, ::Tuple{Symbol,Symbol}, ::Tuple{Nothing,Nothing}) at /home/yakir/.julia/packages/GridLayoutBase/Y582M/src/layoutobservables.jl:267
 [74] (::Observables.OnUpdate{GridLayoutBase.var"#115#116"{Observables.Observable{Tuple{Any,Any}},Observables.Observable{Tuple{Union{Nothing, Float32},Union{Nothing, Float32}}},Observables.Observable{GridLayoutBase.AlignMode},Observables.Observable{GridLayoutBase.RectSides{Float32}},Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}},Tuple{Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}},Observables.Observable{Any},Observables.Observable{Tuple{Union{Nothing, Float32},Union{Nothing, Float32}}}}})(::GeometryBasics.HyperRectangle{2,Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [75] setindex!(::Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}, ::GeometryBasics.HyperRectangle{2,Float32}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [76] setindex!(::Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}, ::GeometryBasics.HyperRectangle{2,Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126
 [77] (::GridLayoutBase.var"#7#9"{GridLayout})(::Bool) at /home/yakir/.julia/packages/GridLayoutBase/Y582M/src/gridlayout.jl:54
 [78] #invokelatest#1 at ./essentials.jl:712 [inlined]
 [79] invokelatest at ./essentials.jl:711 [inlined]
 [80] setindex!(::Observables.Observable{Bool}, ::Bool; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:132
 [81] setindex!(::Observables.Observable{Bool}, ::Bool) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126
 [82] (::GridLayoutBase.var"#14#15"{GridLayout})(::Bool) at /home/yakir/.julia/packages/GridLayoutBase/Y582M/src/gridlayout.jl:138
 [83] #invokelatest#1 at ./essentials.jl:712 [inlined]
 [84] invokelatest at ./essentials.jl:711 [inlined]
 [85] setindex!(::Observables.Observable{Bool}, ::Bool; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:132
 [86] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [87] (::GridLayoutBase.var"#10#12"{GridLayoutBase.GridContent{GridLayout,LAxis}})(::GridLayoutBase.RectSides{Float32}) at /home/yakir/.julia/packages/GridLayoutBase/Y582M/src/gridlayout.jl:109
 [88] #invokelatest#1 at ./essentials.jl:712 [inlined]
 [89] invokelatest at ./essentials.jl:711 [inlined]
 [90] setindex!(::Observables.Observable{GridLayoutBase.RectSides{Float32}}, ::GridLayoutBase.RectSides{Float32}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:132
 [91] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [92] #107 at /home/yakir/.julia/packages/GridLayoutBase/Y582M/src/layoutobservables.jl:23 [inlined]
 [93] (::Observables.OnUpdate{GridLayoutBase.var"#107#109"{Observables.Observable{GridLayoutBase.RectSides{Float32}}},Tuple{Observables.Observable{GridLayoutBase.RectSides{Float32}},Observables.Observable{Any}}})(::GridLayoutBase.RectSides{Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [94] setindex!(::Observables.Observable{GridLayoutBase.RectSides{Float32}}, ::GridLayoutBase.RectSides{Float32}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [95] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [96] (::MakieLayout.var"#139#163"{Observables.Observable{GridLayoutBase.RectSides{Float32}},MakieLayout.var"#compute_protrusions#162"{Text{...}}})(::String, ::Vararg{Any,N} where N) at /home/yakir/.julia/packages/MakieLayout/bYPho/src/lobjects/laxis.jl:300
 [97] (::Observables.OnUpdate{MakieLayout.var"#139#163"{Observables.Observable{GridLayoutBase.RectSides{Float32}},MakieLayout.var"#compute_protrusions#162"{Text{...}}},Tuple{Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Float32},Observables.Observable{Float32},Observables.Observable{Any},Observables.Observable{Any}}})(::Float32) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [98] setindex!(::Observables.Observable{Float32}, ::Float32; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [99] setindex!(::Observables.Observable{Float32}, ::Float32) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126
 [100] (::Observables.MapUpdater{MakieLayout.var"#96#116"{Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Text{...}},Float32})(::Bool, ::Vararg{Any,N} where N) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:241
 [101] (::Observables.OnUpdate{Observables.MapUpdater{MakieLayout.var"#96#116"{Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Text{...}},Float32},Tuple{Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Float32},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any}}})(::Float32) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [102] setindex!(::Observables.Observable{Any}, ::Float32; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [103] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [104] (::MakieLayout.var"#79#99"{Observables.Observable{Any}})(::Float32, ::AbstractPlotting.Automatic) at /home/yakir/.julia/packages/MakieLayout/bYPho/src/lineaxis.jl:68
 [105] (::Observables.OnUpdate{MakieLayout.var"#79#99"{Observables.Observable{Any}},Tuple{Observables.Observable{Float32},Observables.Observable{Any}}})(::Float32) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [106] setindex!(::Observables.Observable{Float32}, ::Float32; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [107] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [108] MapUpdater at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:241 [inlined]
 [109] (::Observables.OnUpdate{Observables.MapUpdater{MakieLayout.var"#78#98"{Observables.Observable{Any},Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Annotations{...}},Float32},Tuple{Observables.Observable{Array{Tuple{String,Point{2,Float32}},1}},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Any}}})(::Array{Tuple{String,Point{2,Float32}},1}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [110] setindex!(::Observables.Observable{Array{Tuple{String,Point{2,Float32}},1}}, ::Array{Tuple{String,Point{2,Float32}},1}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [111] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [112] (::MakieLayout.var"#89#109"{Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Observables.Observable{Array{Tuple{String,Point{2,Float32}},1}},Observables.Observable{Float32},Observables.Observable{Array{Point{2,Float32},1}}})(::Array{String,1}, ::Float32, ::Bool) at /home/yakir/.julia/packages/MakieLayout/bYPho/src/lineaxis.jl:199
 [113] (::Observables.OnUpdate{MakieLayout.var"#89#109"{Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Observables.Observable{Array{Tuple{String,Point{2,Float32}},1}},Observables.Observable{Float32},Observables.Observable{Array{Point{2,Float32},1}}},Tuple{Observables.Observable{Array{String,1}},Observables.Observable{Float32},Observables.Observable{Any}}})(::Array{String,1}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [114] setindex!(::Observables.Observable{Array{String,1}}, ::Array{String,1}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [115] setindex!(::Observables.Observable{Array{String,1}}, ::Array{String,1}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126
 [116] (::MakieLayout.var"#86#106"{Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Observables.Observable{Array{Point{2,Float32},1}},Observables.Observable{Array{String,1}}})(::Array{Float64,1}) at /home/yakir/.julia/packages/MakieLayout/bYPho/src/lineaxis.jl:178
 [117] (::Observables.OnUpdate{MakieLayout.var"#86#106"{Observables.Observable{Any},Observables.Observable{Any},Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Observables.Observable{Array{Point{2,Float32},1}},Observables.Observable{Array{String,1}}},Tuple{Observables.Observable{Array{Float64,1}}}})(::Array{Float64,1}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [118] setindex!(::Observables.Observable{Array{Float64,1}}, ::Array{Float64,1}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [119] setindex!(::Observables.Observable{Array{Float64,1}}, ::Array{Float64,1}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126
 [120] (::Observables.MapUpdater{MakieLayout.var"#85#105",Array{Float64,1}})(::Tuple{Float32,Tuple{Float32,Float32},Bool}, ::Vararg{Any,N} where N) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:241
 [121] (::Observables.OnUpdate{Observables.MapUpdater{MakieLayout.var"#85#105",Array{Float64,1}},Tuple{Observables.Observable{Tuple{Float32,Tuple{Float32,Float32},Bool}},Observables.Observable{Any},Observables.Observable{Any}}})(::Tuple{Float32,Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [122] setindex!(::Observables.Observable{Any}, ::Tuple{Float32,Float32}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [123] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [124] (::Observables.var"#3#4"{Any,Observables.Observable{Any}})(::Tuple{Float32,Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:47
 [125] #invokelatest#1 at ./essentials.jl:712 [inlined]
 [126] invokelatest at ./essentials.jl:711 [inlined]
 [127] setindex!(::Observables.Observable{Tuple{Float32,Float32}}, ::Tuple{Float32,Float32}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:132
 [128] setindex! at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126 [inlined]
 [129] MapUpdater at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:241 [inlined]
 [130] (::Observables.OnUpdate{Observables.MapUpdater{typeof(MakieLayout.ylimits),Tuple{Float32,Float32}},Tuple{Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}}})(::GeometryBasics.HyperRectangle{2,Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [131] setindex!(::Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}, ::GeometryBasics.HyperRectangle{2,Float32}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [132] setindex!(::Observables.Observable{GeometryBasics.HyperRectangle{2,Float32}}, ::GeometryBasics.HyperRectangle{2,Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126
 [133] adjustlimits!(::LAxis) at /home/yakir/.julia/packages/MakieLayout/bYPho/src/lobjects/laxis.jl:603
 [134] #140 at /home/yakir/.julia/packages/MakieLayout/bYPho/src/lobjects/laxis.jl:329 [inlined]
 ... (the last 83 lines are repeated 331 more times)
 [27608] (::Observables.OnUpdate{MakieLayout.var"#140#164"{LAxis},Tuple{Observables.Observable{GeometryBasics.HyperRectangle{2,Int64}},Observables.Observable{Any}}})(::GeometryBasics.HyperRectangle{2,Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:218
 [27609] setindex!(::Observables.Observable{Any}, ::GeometryBasics.HyperRectangle{2,Float32}; notify::Observables.var"#6#8") at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:130
 [27610] setindex!(::Observables.Observable{Any}, ::GeometryBasics.HyperRectangle{2,Float32}) at /home/yakir/.julia/packages/Observables/0wrF6/src/Observables.jl:126
 [27611] autolimits!(::LAxis) at /home/yakir/.julia/packages/MakieLayout/bYPho/src/lobjects/laxis.jl:549
 [27612] plot!(::LAxis, ::Type{Poly{...}}, ::Attributes, ::Array{Point{2,Float64},1}; kw_attributes::Base.Iterators.Pairs{Union{},Union{},Tuple{},NamedTuple{(),Tuple{}}}) at /home/yakir/.julia/packages/MakieLayout/bYPho/src/lobjects/laxis.jl:393
 [27613] plot! at /home/yakir/.julia/packages/MakieLayout/bYPho/src/lobjects/laxis.jl:391 [inlined]
 [27614] #poly!#307 at /home/yakir/.julia/packages/AbstractPlotting/WwXMh/src/recipes.jl:24 [inlined]
 [27615] top-level scope at /home/yakir/coffeebeetlearticle/displacement/main.jl:135
 [27616] include(::String) at ./client.jl:439
in expression starting at /home/yakir/coffeebeetlearticle/displacement/main.jl:91
