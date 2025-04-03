using DynamicalSystems
include("ComplexDynamicalSystems.jl")
using .ComplexDynamicalSystems


struct EscapeTimeConfiguration
    parameter_grid::ComplexGridConfiguration
    max_dist::Float64
    max_steps::Int
end

function escapes(system::PolyHolomorphicDS, configuration::EscapeTimeConfiguration)
    for t in 1:configuration.max_steps
        step!(system.state)
        if abs(fpf(current_state(system.state))) > configuration.max_dist
            return true
        end
    end
    return false
end

function BuildMandelbrotSet(system::Any, configuration::EscapeTimeConfiguration)
    parameter_grid = ComplexGrid(configuration.parameter_grid, false)
    for i in 1:size(parameter_grid.elems, 1)
        for j in 1:size(parameter_grid.elems, 2)
            if !escapes(system(parameter_grid.elems[i,j]), configuration)  
                parameter_grid.marks[i, j] = true
            end
        end
    end
    return parameter_grid
end