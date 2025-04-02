

module ComplexDynamicalSystems

using StaticArrays 
using DynamicalSystems
using Polynomials

export tpf, fpf
export PolyHolomorphicDS
export ComplexVariable, ComplexPolynomial, inverse
export InversionBranchMethod, full, random
export ComplexGrid, ComplexGridConfiguration, gridmk, snap_to_grid 

include("ComplexDynamicalSystems/utills.jl")
include("ComplexDynamicalSystems/complex-polynomials.jl")
include("ComplexDynamicalSystems/complex-systems.jl")


end
