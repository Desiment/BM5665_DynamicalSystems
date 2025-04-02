function tpf(z::ComplexF64)
    return z.re, z.im
end
function fpf(x::Float64, y::Float64)
    return x + y*1im
end
function fpf(coords::SVector{2, Float64})
    return fpf(coords[1],coords[2])
end

function _check_grid(ld::ComplexF64, ru::ComplexF64)
    return !((ld.re < ru.re) && (ld.im < ru.im))
end

struct ComplexGridConfiguration
    resolution::Int
    ld::ComplexF64
    lu::ComplexF64
    ru::ComplexF64
    rd::ComplexF64
    ComplexGridConfiguration(r::Int, ld::ComplexF64, ru::ComplexF64) = 
        _check_grid(ld, ru) ? error("incorrect grid") : new(r, ld, ld.re + ru.im*1im,
                                                               ru, ru.re + ld.im*1im)
end

mutable struct ComplexGrid{T}
    configuration::ComplexGridConfiguration
    marks::Matrix{T}
    elems::Matrix{ComplexF64}
    ComplexGrid(cfg::ComplexGridConfiguration, val) = new{typeof(val)}(cfg, fill(val, size(gridmk(cfg))), gridmk(cfg)) 
end


function inside_grid_border(grid::ComplexGridConfiguration, pt::ComplexF64)
    return ((grid.ld.re <= pt.re) && (pt.re <= grid.ru.re) &&
            (grid.ld.im <= pt.im) && (pt.im <= grid.ru.im))
end

function gridmk(g::ComplexGridConfiguration)
    return [i + j*1.0im  for i in g.ld.re:(10.0^(-g.resolution)):g.ru.re, j in g.ld.im:(10.0^(-g.resolution)):g.ru.im]
end 

function snap_to_grid(g::ComplexGridConfiguration, values::Vector{ComplexF64})
    return filter(x->inside_grid_border(g,x),
        [round(values[i], digits=g.resolution)  for i in 1:length(values)])
end

Base.:^(f::Union{Type,Function}, n::Integer) = n <= 1 ? f : f ∘ f^(n-1) # little hack

function bitmap_to_file(grid::ComplexGrid{Bool}, border::Int, c_inverse::Bool, filename::String)
    border_color = c_inverse ? 0 : 1
    open(filename, "w") do f
        for i in 1:border
            for j in 1:(2 * border + size(grid.marks,2))
                print(f, border_color, " ")
            end
            println(f, " ")
        end
        for i in 1:size(grid.marks,1)
            for j in 1:border
                print(f, border_color, " ")
            end
            for j in  1:size(grid.marks,2)
                print(f, grid.marks[i,j] ? (1-border_color) : border_color, " ")
            end
            for j in 1:border
                print(f, border_color, " ")
            end
            println(f, " ")
        end
        for i in 1:border
            for j in 1:(2 * border + size(grid.marks,2))
                print(f, border_color, " ")
            end
            println(f, " ")
        end
    end   
end