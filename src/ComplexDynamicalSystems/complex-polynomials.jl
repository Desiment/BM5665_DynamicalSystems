
ComplexPolynomial = Polynomial{ComplexF64}
(poly::ComplexPolynomial)(points::Vector{ComplexF64}) = [poly(points[i]) for i in 1:length(points)]


function ComplexVariable()
    return  ComplexPolynomial([0.0, 1.0])
end

@enum InversionBranchMethod begin
    full   = 1
    random = 2
end

function inverse(poly::ComplexPolynomial, point::ComplexF64, method::InversionBranchMethod)
    poly_roots = []
    if degree(poly) == 2
        D = (poly[1])^2 - 4*(poly[0] - point)*poly[2]
        poly_roots = [((poly[1]) + sqrt(D)) / (2 * poly[2]), ((poly[1]) - sqrt(D)) / (2 * poly[2])]
    else
        poly_roots = roots(poly - ComplexPolynomial([point]))
    end
    
    if method == full
        return poly_roots
    elseif method == random
        return rand(poly_roots)
    end
end

function inverse(poly::ComplexPolynomial, points::Vector{ComplexF64}, method::InversionBranchMethod)
    return reduce(vcat, [inverse(poly, points[i], method) for i in 1:length(points)])
end
