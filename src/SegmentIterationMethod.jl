Point = Tuple{Float64, Float64};
PolygonalChain = Vector{Point};

struct SegmentIterrationConfiguration
    domain::Pair{Point}
    niters::Int
    fragm_dist::Real
    chain_base::Vector{Point}
end

function distance(a::Point, b::Point)
    return ((a[1]-b[1])^2 + (a[2]-b[2])^2)^0.5
end

function inside(a::Point, domain::Pair{Point})
    return (domain[1][1] <= a[1]) & (a[1] <= domain[2][1]) &
           (domain[1][2] <= a[2]) & (a[2] <= domain[2][2])
end

function neatuniq(v)
    v1 = Vector{eltype(v)}()
    if length(v) == 0
        return v1
    end
    last = v[1]
    push!(v1,last)
    for e in v
        if e != last
            last = e
            push!(v1,last)
        end
    end
    return v1
end

function fragmentize(chain::PolygonalChain, system::Any, configuration::SegmentIterrationConfiguration)
    fragmentized = Vector{Point}()
    for i in 1:(length(chain)-1)
        if !(inside(system(chain[i]), configuration.domain) & inside(system(chain[i+1]), configuration.domain))
            continue
        end
        if distance(system(chain[i]), system(chain[i+1])) < configuration.fragm_dist
            push!(fragmentized, chain[i])
            push!(fragmentized, chain[i+1])
        else
            midpoint = (chain[i] .+ chain[i+1]) ./ 2
            append!(fragmentized, fragmentize([chain[i], midpoint], system, configuration))
            append!(fragmentized, fragmentize([midpoint, chain[i+1]], system, configuration))
        end
    end
	return neatuniq(fragmentized)
end

function SegmentIterration(system::Any, configuration::SegmentIterrationConfiguration)
    chain = configuration.chain_base
    for i in 1:configuration.niters
	    chain = fragmentize(chain, system, configuration)
        chain = system.(chain)
    end
    return chain
end
