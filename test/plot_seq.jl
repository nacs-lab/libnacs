#!/usr/bin/julia

using PyPlot

function add_point(dict, chn_name, chn_id, t, v)
    key = chn_name, chn_id
    if key in keys(dict)
        ts, vs = dict[key]
        v_prev = vs[end]
        push!(ts, t)
        push!(ts, t)
        push!(vs, v_prev)
        push!(vs, v)
    else
        ts, vs = dict[key] = Float64[], Float64[]
        push!(ts, t)
        push!(vs, v)
    end
    return
end

function load_seq_txt(fname,
                      dict=Dict{Tuple{String,Int},
                                Tuple{Vector{Float64},Vector{Float64}}}())
    open(fname) do fh
        for line in eachline(fh)
            m = match(r"^([0-9]*), ([^()]*)\(([^()]*)\) = (.*)",
                      line)::RegexMatch
            t = parse(Int, m[1]) * 10e-9
            chn_name = m[2]
            if chn_name == "TTL"
                @assert m[3] == "all"
                vttl = parse(UInt32, m[4])
                for i in 0:31
                    add_point(dict, "TTL", i, t, ((UInt32(1) << i) & vttl) != 0)
                end
                # @printf "%.9f, TTL, all, 0x%08x\n" t vttl
            else
                chn_id = parse(Int, m[3])
                vf64 = parse(Float64, m[4])
                add_point(dict, String(chn_name), chn_id, t, vf64)
                # @printf "%.9f, %s, %d, %.8e\n" t chn_name chn_id vf64
            end
        end
    end
    dict
end

function filter_data(dict, chn::String)
    dict2 = Dict{Tuple{String,Int},Tuple{Vector{Float64},Vector{Float64}}}()
    for (k, v) in dict
        if k[1] == chn
            dict2[k] = v
        end
    end
    dict2
end

function filter_data(dict, chn::Tuple)
    dict2 = Dict{Tuple{String,Int},Tuple{Vector{Float64},Vector{Float64}}}()
    for (k, v) in dict
        if k == chn
            dict2[k] = v
        end
    end
    dict2
end

function plot_channel(dict, chn)
    dict2 = filter_data(dict, chn)
    for (k, v) in dict2
        chn_name = "$(k[1])($(k[2]))"
        plot(v[1], v[2], label=chn_name)
    end
end

# display(load_seq_txt(ARGS[1]))
