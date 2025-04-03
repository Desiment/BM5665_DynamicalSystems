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

# ╔═╡ 6de23c3c-8be5-4e8b-a24d-826f63e0e5d3
# ╠═╡ show_logs = false
begin
	import Pkg; Pkg.activate("..") #use project's manifest
	#imports
	using PlutoUI, HypertextLiteral, Plots
	using ImageShow, ImageIO, Colors
	using PlutoImageCoordinatePicker
	import PlutoUI: combine 
end

# ╔═╡ f2b458a6-100c-11f0-0037-131a512b6134
md"# Задача VI. Модель логистической решетки"

# ╔═╡ f2b4591e-100c-11f0-139d-131cbcc569f4
md"Модель логистической решетки реализована ниже (считаем что стороны квадрата склеены в тор):"

# ╔═╡ f2b4598c-100c-11f0-34a4-491c00863324
md"Возьмем квадрат 30x30, в качестве параметров возьмем $\mathrm{GN} = 900, \mathrm{CP} = 500$"

# ╔═╡ dd005166-2b85-4ee5-b1ce-78da83e291de
@bind llmodel_configuration confirm(
	combine() do Child
	@htl("""
		<h5>Параметры метода</h5>
		<ul>
		$([
			@htl("<li> Размер решетки, n = $(Child("n", NumberField(10:100; default=30)))"),
			@htl("<li>Параметры модели. GN = $(Child("GN", NumberField(0:2000; default=900))), CP =  $(Child("CP", NumberField(0:2000; default=500)))"),
			@htl("<li>Число итераций: $(Child("t", NumberField(0:200; default=100)))")
		])
		</ul>
		""")
	end	
	)

# ╔═╡ f2b459a0-100c-11f0-2420-7ff021cefd4e
begin
	u0 = rand(llmodel_configuration.n, llmodel_configuration.n)
	p0 = [4 * llmodel_configuration.GN / 1000, llmodel_configuration.CP / 1000];
	states = [u0]
end

# ╔═╡ f2b459aa-100c-11f0-3c53-bf58f42b0436
md"Произведем моделирование"

# ╔═╡ 5f301fcd-f815-4ee6-ab2e-52e746957859
md"""
# Прочее
"""

# ╔═╡ f2b459c6-100c-11f0-04fa-cd26e97c83f6
function trim(x::Real, a::Real, b::Real)
    return x < a ? a : (x > b ? b : x)
end

# ╔═╡ f2b4596e-100c-11f0-200f-3987d236356d
function logistic_lattice_step(u, p)
    averaged = (circshift(u, (0, -1)) .+ circshift(u, (-1, 0)) .+ circshift(u, (1, 0)) .+ circshift(u, (0, 1))) ./ 4
    state = p[1] .* u .* (1 .- u) .+ p[2] .* (averaged - u)
    return trim.(state, 0, 1)
end

# ╔═╡ f2b459be-100c-11f0-25d9-cd1d05e53b26
anim = @animate for i in 1:llmodel_configuration.t
	    global states
	    push!(states, logistic_lattice_step(states[end], p0))
	    heatmap(states[end], clim=(0,1))
	    title!(("Logistic lattice, step "*string(i)))
end

# ╔═╡ 2e626f1c-1b0b-41ae-99da-fb845e884f4f
gif(anim, "results/problem-06-logistic-lattice.gif", fps = 1)

# ╔═╡ Cell order:
# ╟─6de23c3c-8be5-4e8b-a24d-826f63e0e5d3
# ╟─f2b458a6-100c-11f0-0037-131a512b6134
# ╟─f2b4591e-100c-11f0-139d-131cbcc569f4
# ╟─f2b4598c-100c-11f0-34a4-491c00863324
# ╟─dd005166-2b85-4ee5-b1ce-78da83e291de
# ╠═f2b4596e-100c-11f0-200f-3987d236356d
# ╠═f2b459a0-100c-11f0-2420-7ff021cefd4e
# ╟─f2b459aa-100c-11f0-3c53-bf58f42b0436
# ╟─f2b459be-100c-11f0-25d9-cd1d05e53b26
# ╟─2e626f1c-1b0b-41ae-99da-fb845e884f4f
# ╟─5f301fcd-f815-4ee6-ab2e-52e746957859
# ╠═f2b459c6-100c-11f0-04fa-cd26e97c83f6
