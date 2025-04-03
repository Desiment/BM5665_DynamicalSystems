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

# ╔═╡ 32c5fc17-9d58-4858-8ebd-e5396099a5e1
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
	include("../src/EscapeTimeMethod.jl") # method implementation
	using .ComplexDynamicalSystems
end

# ╔═╡ eef5c592-100c-11f0-2339-41b1ab2fbcd8
md"# Задача IV. Кубическое множество Мандельброта"

# ╔═╡ eefa3e42-100c-11f0-3969-5366f782635a
md"""
Построить множество Мандельброта для отображения $f(z) = z^3 + c$.

### Решение
Для построения используется метод выхода за границу области (escape time algorithm). 

* Вспомогательные методы реализованы в `src/ComplexDynamicalSystems.jl`. 
* Метод реализован в файле `src/EscapeTimeMethod.jl`
"""

# ╔═╡ eefaa5d0-100c-11f0-2533-a5e57fa26972
md"Определим семейство отображений как функцию от параметра:"

# ╔═╡ eefaa634-100c-11f0-00a7-e596df5557db
function CubicFamily(c::ComplexF64)
    z = ComplexVariable()
    return PolyHolomorphicDS(c + z^3)
end;

# ╔═╡ eefaa666-100c-11f0-0870-61a1d88b0bdf
md"""
Конфигурация метода:
- Число итераций для каждого параметра: $n$
- Область, при выходе за границу которой, параметр отбрасывается из множества: $|z| \leq R$
- Рассматриваются параметры лежащие в квадратной сетке : $\max(|\mathrm{Re} z|, |\mathrm{Im} z|) \leq d$ с шагом $h$ по вещественной и мнимым осям"""

# ╔═╡ ac494e9f-b9bb-4394-94c0-3498f62b5fdd
@bind etm_parameters confirm(
	combine() do Child
	@htl("""
		<h5>Параметры метода</h5>
		<ul>
		$([
			@htl("<li> Радиус ухода, R = $(Child("R", NumberField(0:50; default=10)))"),
			@htl("<li> Сетка параметров. Параметр границы сетки, d = $(Child("d", NumberField(0:0.1:2; default=1))), шаг h =  10^-$(Child("hexp", NumberField(1:4; default=3)))"),
			@htl("<li>Число итераций: $(Child("n", NumberField(0:500; default=100)))")
		])
		</ul>
		""")
	end	
	)

# ╔═╡ c15ea9c3-ed1d-4c0b-bb65-c820293581f7
@bind output_filename confirm(@htl("Choose output filename $(TextField())"))

# ╔═╡ 00a6ce1f-3d2a-49be-a149-7d6c26c07e12
md"""
# Прочее
"""

# ╔═╡ 8ad16e4a-2180-422b-9a21-38429165b311
function etm_configuration_from_ui(method_params)
	return EscapeTimeConfiguration(
			ComplexGridConfiguration(method_params.hexp, 
							-method_params.d - method_params.d*1im,
							 method_params.d + method_params.d*1im), # parameter_grid 
			method_params.R, # max_distance
		    method_params.n # max_itters
			)
end

# ╔═╡ eefaa682-100c-11f0-3c37-19e41d48cef8
mdb = BuildMandelbrotSet(CubicFamily, etm_configuration_from_ui(etm_parameters));

# ╔═╡ eefaa68e-100c-11f0-1080-93cf0acc1fd1
begin
	mdb_set=scatter(real(mdb.elems[mdb.marks]), imag(mdb.elems[mdb.marks]), m=0, ms=0, title="Mandelbrot set for f(z) = z³ + c", legend=false, aspect_ratio=:equal);
	xlabel!("Re")
	ylabel!("Im")
end

# ╔═╡ 5404378c-b517-4816-95f5-8311d10e7cdb
savefig(mdb_set, if output_filename != "" "results/"*output_filename else "results/problem-04" end)

# ╔═╡ eefaa698-100c-11f0-23e6-0f91da5981cb
md"![](results/T4_cubic_mandelbrot_set.png)"

# ╔═╡ Cell order:
# ╟─32c5fc17-9d58-4858-8ebd-e5396099a5e1
# ╟─eef5c592-100c-11f0-2339-41b1ab2fbcd8
# ╟─eefa3e42-100c-11f0-3969-5366f782635a
# ╟─eefaa5d0-100c-11f0-2533-a5e57fa26972
# ╠═eefaa634-100c-11f0-00a7-e596df5557db
# ╟─eefaa666-100c-11f0-0870-61a1d88b0bdf
# ╟─ac494e9f-b9bb-4394-94c0-3498f62b5fdd
# ╠═eefaa682-100c-11f0-3c37-19e41d48cef8
# ╠═eefaa68e-100c-11f0-1080-93cf0acc1fd1
# ╠═c15ea9c3-ed1d-4c0b-bb65-c820293581f7
# ╟─5404378c-b517-4816-95f5-8311d10e7cdb
# ╟─00a6ce1f-3d2a-49be-a149-7d6c26c07e12
# ╟─8ad16e4a-2180-422b-9a21-38429165b311
# ╟─eefaa698-100c-11f0-23e6-0f91da5981cb
