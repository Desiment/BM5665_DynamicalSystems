include("ComplexDynamicalSystems.jl")
using .ComplexDynamicalSystems

struct InvIterMethodConfiguration
    niters::Int
    branch::InversionBranchMethod
    grid::ComplexGridConfiguration
end

function BuildJuliaSet(system::PolyHolomorphicDS, 
                       config::InvIterMethodConfiguration)
    return BuildJuliaSet(system, config, reduce(vcat, gridmk(config.grid)))
end

function BuildJuliaSet(system::PolyHolomorphicDS, 
                       config::InvIterMethodConfiguration,
                       init::Vector{ComplexF64})
    base_grid = ComplexGrid(config.grid, false)
    points = init
    for _ in 1:config.niters
        points = inverse(system.map, points, config.branch)
        points = snap_to_grid(config.grid, points)
        points = unique(points)
    end
    idxs = map(identity, filter(x -> !isnothing(x), indexin(points, base_grid.elems)))
    base_grid.marks[idxs] .= true
    return base_grid
end
