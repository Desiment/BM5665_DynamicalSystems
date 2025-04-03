### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a7a20bc8-fef1-4e08-85e2-842fe4930c00
# ╠═╡ show_logs = false
begin
	import Pkg; Pkg.activate("..") #use project's manifest
	#imports
	using PlutoUI, HypertextLiteral, Plots
	using ImageShow, ImageIO, Colors
	using PlutoImageCoordinatePicker
	import PlutoUI: combine 
	# sources
	include("../src/ComplexDynamicalSystems.jl") # basis utills
	include("../src/FindPeriodicOrbits.jl") # method implementation
	using .ComplexDynamicalSystems
	using Polynomials
end

# ╔═╡ f100913c-100c-11f0-010e-c7f6c3d1cf75
md"# Задача V. Периодические орбиты функции Ньютона."

# ╔═╡ f1009218-100c-11f0-0a91-c7a989b64fdb
md"""
Найти приближения к периодическим орбитам функции Ньютона для многочлена $f(z) = z^3 + (A - 1)z - A$. Значения параметра:  
*  $A = 0.310 + 1.620i$ (период 2);
*  $A = 0.275 + 1.650i$ (период 4).
"""

# ╔═╡ f10092a4-100c-11f0-3d5c-6db969c77f5c
md"
### Решение
Ниже строятся приближения к множествам точек таких что $N_f^p(z) = z$ (т.е. например для $p = 4$ в эти множества попадают также неподвижные точки и орбиты периода 2), где 

$N_f(z) = z - f(z)/f'(z)$

функция Ньютона для $f(z)$, $N^p(z) = (N \circ \ldots \circ N)(z)$, $p$ раз"

# ╔═╡ f10092c4-100c-11f0-0c01-7f53516912ae
md"""##### Первый случай, орбиты периода 2

Объявим необходимые сущности:
"""

# ╔═╡ f10092e0-100c-11f0-0959-17ec3246cde2
begin
	z = ComplexVariable()
	A1 = 0.310 + 1.620im
	f1 = -A1 + (A1 - 1) * z + z^3;
end

# ╔═╡ f10092fe-100c-11f0-1b58-119173242eea
md"Определим функцию Ньютона:"

# ╔═╡ f1009312-100c-11f0-0b8c-cdd723b4bb13
newton_c1(x) = (x - f1(x) / derivative(f1)(x));

# ╔═╡ f1009324-100c-11f0-0ba4-8382d3cc1d3d
md"Искать будем на сетке  $\max(|\mathrm{Re} z|, |\mathrm{Im} z|) \leq 1.5$ с шагом $0.001$ по вещественной и мнимым осям"

# ╔═╡ f100933a-100c-11f0-20ed-ef08046c7c53
grid_configuration = ComplexGridConfiguration(3, -1.5 - 1.5im, 1.5 + 1.5im);

# ╔═╡ f1009356-100c-11f0-3e70-abcc39f6d9a1
md"Ищем периодические орбиты порядка 2, проверяя для каждой точки на сетке что $|N^2(z) - z| < \varepsilon$, где $\varepsilon = 0.1, 0.075, 0.05, 0.025, 0.01, 0.001$"

# ╔═╡ f1009380-100c-11f0-358d-2b0fb01d9616
begin
	period_c1 = 2
	accuracy_c1 = [0.1, 0.075, 0.05, 0.025, 0.01, 0.001]
	approximations_c1 = Vector{ComplexGrid{Bool}}()
	for ε in accuracy_c1
	    search_grid = ComplexGrid(grid_configuration, false)
	    push!(approximations_c1, find_approx_periodic_orbits(newton_c1, period_c1, search_grid, ε))
	end
end

# ╔═╡ f1009394-100c-11f0-02c6-37d328c4a7ef
md"Построим полученные множества:"

# ╔═╡ f10093da-100c-11f0-09b3-2338940b4030
begin
l_c1 = @layout [a{0.04h}; grid(2,3)]
p1_c1 = scatter(real(approximations_c1[1].elems[approximations_c1[1].marks]), 
             imag(approximations_c1[1].elems[approximations_c1[1].marks]),
             m=0, ms=0, title="ε = 0.1", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
p2_c1 = scatter(real(approximations_c1[2].elems[approximations_c1[2].marks]), 
             imag(approximations_c1[2].elems[approximations_c1[2].marks]),
             m=0, ms=0, title="ε = 0.075", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
p3_c1 = scatter(real(approximations_c1[3].elems[approximations_c1[3].marks]), 
             imag(approximations_c1[3].elems[approximations_c1[3].marks]),
             m=0, ms=0, title="ε = 0.050", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
p4_c1 = scatter(real(approximations_c1[4].elems[approximations_c1[4].marks]), 
             imag(approximations_c1[4].elems[approximations_c1[4].marks]),
             m=0, ms=0, title="ε = 0.025", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
p5_c1 = scatter(real(approximations_c1[5].elems[approximations_c1[5].marks]), 
             imag(approximations_c1[5].elems[approximations_c1[5].marks]),
             m=0, ms=0, title="ε = 0.010", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
p6_c1 = scatter(real(approximations_c1[6].elems[approximations_c1[6].marks]), 
             imag(approximations_c1[6].elems[approximations_c1[6].marks]),
             title="ε = 0.001", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
title_c1 = plot(title = "Approximation of periodic orbits of order 2 for\n Newton's function of f(z) = z³ + (A - 1)z - A, A = $A1", grid = false, showaxis = false, bottom_margin = -50Plots.px)
c1_plot = plot(title_c1, p1_c1, p2_c1, p3_c1, p4_c1, p5_c1, p6_c1, layout = l_c1);
plot!(size=(1600,1000))
end
#savefig("results/T5_period2_A_0.310_1.620.png");

# ╔═╡ 0e513fb9-1461-42cb-9b38-445a90b6829e
@bind output_filename_c1 confirm(@htl("Choose output filename $(TextField())"))

# ╔═╡ de754bf5-43c9-4121-8190-8ade1a284e10
savefig(c1_plot, "results/"*output_filename_c1)

# ╔═╡ f10093ee-100c-11f0-26bf-a91b9bb3e50e
md"![](results/T5_period2_A_0.310_1.620.png)"

# ╔═╡ f1009402-100c-11f0-358a-2b6cfdcefffd
md"##### Второй случай, орбиты порядка 4"

# ╔═╡ f1009420-100c-11f0-01a5-f302530678e0
begin
	A2 = 0.275 + 1.650im
	f2 = -A2 + (A2 - 1) * z + z^3
	newton_c2(x) = (x - f2(x) / derivative(f2)(x))
end

# ╔═╡ f1009434-100c-11f0-06cb-878c366ffefa
md"Искать будем на той же сетке  $\max(|\mathrm{Re} z|, |\mathrm{Im} z|) \leq 1.5$ с шагом $0.001$ по вещественной и мнимым осям с теми же погрешностями"

# ╔═╡ f1009470-100c-11f0-356c-d7aaf86aeccc
begin
	period_c2 = 4
	accuracy_c2 = [0.1, 0.075, 0.05, 0.025, 0.01, 0.001]
	approximations_c2 = Vector{ComplexGrid{Bool}}()
	for ε in accuracy_c2
	    search_grid = ComplexGrid(grid_configuration, false)
	    push!(approximations_c2, find_approx_periodic_orbits(newton_c2, period_c2, search_grid, ε))
	end
end

# ╔═╡ f10094a2-100c-11f0-22c1-e7b6428828cb
begin
l_c2 = @layout [a{0.04h}; grid(2,3)]
p1_c2 = scatter(real(approximations_c2[1].elems[approximations_c2[1].marks]), 
             imag(approximations_c2[1].elems[approximations_c2[1].marks]),
             m=0, ms=0, title="ε = 0.1", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
p2_c2 = scatter(real(approximations_c2[2].elems[approximations_c2[2].marks]), 
             imag(approximations_c2[2].elems[approximations_c2[2].marks]),
             m=0, ms=0, title="ε = 0.075", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
p3_c2 = scatter(real(approximations_c2[3].elems[approximations_c2[3].marks]), 
             imag(approximations_c2[3].elems[approximations_c2[3].marks]),
             m=0, ms=0, title="ε = 0.050", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
p4_c2 = scatter(real(approximations_c2[4].elems[approximations_c2[4].marks]), 
             imag(approximations_c2[4].elems[approximations_c2[4].marks]),
             m=0, ms=0, title="ε = 0.025", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
p5_c2 = scatter(real(approximations_c2[5].elems[approximations_c2[5].marks]), 
             imag(approximations_c2[5].elems[approximations_c2[5].marks]),
             m=0, ms=0, title="ε = 0.010", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
p6_c2 = scatter(real(approximations_c2[6].elems[approximations_c2[6].marks]), 
             imag(approximations_c2[6].elems[approximations_c2[6].marks]),
             title="ε = 0.001", xlabel = "Re", ylabel = "Im", legend=false, bottom_margin=50*Plots.px);
title_c2 = plot(title = "Approximation of periodic orbits of order 4 for\n Newton's function of f(z) = z³ + (A - 1)z - A, A = $A2", grid = false, showaxis = false, bottom_margin = -50Plots.px)
c2_plot = plot(title_c2, p1_c2, p2_c2, p3_c2, p4_c2, p5_c2, p6_c2, layout = l_c2);
plot!(size=(1600,1200))
end

# ╔═╡ 60eefe01-44ac-41d4-87ca-c66c08581177
@bind output_filename_c2 confirm(@htl("Choose output filename $(TextField())"))

# ╔═╡ 8ffb4cb4-bda8-44a7-a16b-d7a30701e7a0
savefig(c2_plot, "results/"*output_filename_c2)

# ╔═╡ f10094b6-100c-11f0-1188-9fdf3e028986
md"![](results/T5_period2_A_0.275_1.650.png)"

# ╔═╡ Cell order:
# ╟─a7a20bc8-fef1-4e08-85e2-842fe4930c00
# ╟─f100913c-100c-11f0-010e-c7f6c3d1cf75
# ╟─f1009218-100c-11f0-0a91-c7a989b64fdb
# ╟─f10092a4-100c-11f0-3d5c-6db969c77f5c
# ╟─f10092c4-100c-11f0-0c01-7f53516912ae
# ╠═f10092e0-100c-11f0-0959-17ec3246cde2
# ╟─f10092fe-100c-11f0-1b58-119173242eea
# ╠═f1009312-100c-11f0-0b8c-cdd723b4bb13
# ╟─f1009324-100c-11f0-0ba4-8382d3cc1d3d
# ╠═f100933a-100c-11f0-20ed-ef08046c7c53
# ╟─f1009356-100c-11f0-3e70-abcc39f6d9a1
# ╠═f1009380-100c-11f0-358d-2b0fb01d9616
# ╟─f1009394-100c-11f0-02c6-37d328c4a7ef
# ╟─f10093da-100c-11f0-09b3-2338940b4030
# ╠═0e513fb9-1461-42cb-9b38-445a90b6829e
# ╟─de754bf5-43c9-4121-8190-8ade1a284e10
# ╟─f10093ee-100c-11f0-26bf-a91b9bb3e50e
# ╟─f1009402-100c-11f0-358a-2b6cfdcefffd
# ╠═f1009420-100c-11f0-01a5-f302530678e0
# ╟─f1009434-100c-11f0-06cb-878c366ffefa
# ╠═f1009470-100c-11f0-356c-d7aaf86aeccc
# ╟─f10094a2-100c-11f0-22c1-e7b6428828cb
# ╟─60eefe01-44ac-41d4-87ca-c66c08581177
# ╟─8ffb4cb4-bda8-44a7-a16b-d7a30701e7a0
# ╟─f10094b6-100c-11f0-1188-9fdf3e028986
