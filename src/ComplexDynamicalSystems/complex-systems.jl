function PolynomialMap(u::SVector{2, Float64}, p::Vector{ComplexF64}, t::Int64)
    return SVector(tpf(ComplexPolynomial(p)(u[1] + u[2] * 1im)))
end

mutable struct PolyHolomorphicDS
    map::ComplexPolynomial
    state::DeterministicIteratedMap
    
    PolyHolomorphicDS(p::Vector{ComplexF64}) = new(ComplexPolynomial(p),
                            DeterministicIteratedMap(PolynomialMap, SVector(Float64(0.0),Float64(0.0)), p))
    PolyHolomorphicDS(p::Vector{ComplexF64}, u::ComplexF64) = new(ComplexPolynomial(p),
                            DeterministicIteratedMap(PolynomialMap, SVector(tpf(u)),  coeffs(p)))
    PolyHolomorphicDS(p::ComplexPolynomial) = new(p,
                            DeterministicIteratedMap(PolynomialMap, SVector(Float64(0.0),Float64(0.0)), coeffs(p)))
    PolyHolomorphicDS(p::ComplexPolynomial, u::ComplexF64) = new(p,
                            DeterministicIteratedMap(PolynomialMap, SVector(tpf(u)), p.coeffs))
end

