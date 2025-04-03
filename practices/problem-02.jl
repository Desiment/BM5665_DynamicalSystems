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

# ╔═╡ 0142c4dc-0fdf-11f0-3403-77788dee9947
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
	include("../src/InverseIterationMethod.jl") # method implementation
	using .ComplexDynamicalSystems
end

# ╔═╡ 013db4c4-0fdf-11f0-0aaf-c92fc9d78215
md"# Задача II. Метод обратных итераций"

# ╔═╡ dd7de283-d8a8-47ed-96a1-98ee4d715abd
md"""
Построить методом обратных итераций приближения к инвариантным множествам отображения $z_{n+1} = z_n^2 + c$, где $z, c \in \mathbb C$, $c$ - параметр.
"""

# ╔═╡ 0142570c-0fdf-11f0-2ec7-1378fe6665f6
md"""
### Решение
* Вспомогательные методы реализованы в `src/ComplexDynamicalSystems.jl`.  
* Метод реализован в файле `src/InverseIterationMethod`
"""

# ╔═╡ 0142c5ac-0fdf-11f0-11ec-6919b42c22a9
md"""
Конфигурация метода:
- Число обратных итераций: $n$
- При обращении рассматриваются обе ветви квадратного корня
- Область в которой рассматривается множество: $\max(|\mathrm{Re} z|, |\mathrm{Im} z|) \leq d$
- Одновременно обрабатываются обратные итерации точек квадратной сетки $\max(|\mathrm{Re} z|, |\mathrm{Im} z|) \leq \ell$ с шагом $h = 0.1$ по вещественной и мнимым осям

"""

# ╔═╡ 9cd321b6-f105-434c-b0ae-5fc44c6a6988
@bind iim_parameters confirm(
	combine() do Child
	@htl("""
		<h5>Параметры метода</h5>
		<ul>
		$([
			@htl("<li> Параметр границы области, d = $(Child("d", NumberField(0:0.1:2; default=1.5)))"),
			@htl("<li>Параметры начальной сетки точек. Параметр границы сетки, l = $(Child("l", NumberField(0:0.1:2; default=1))), шаг h =  $(Child("h", NumberField(0:0.05:1; default=0.1)))"),
			@htl("<li>Число итераций: $(Child("n", NumberField(0:50; default=40)))")
		])
		</ul>
		""")
	end	
	)

# ╔═╡ 0142c55c-0fdf-11f0-1369-d1e8b81243e6
md"
### Применение процедуры
Для примера построим множество Жюлиа для динамической системы, порождаемой отображением $f(z) = c + z^2, c = 0.285$ (в случае интерактивного запуска параметр $c$ можно выбрать с помощью слайдера ниже). Объявим эту систему"

# ╔═╡ bfc236f3-68ed-4da4-8fcf-8cea6363235d
md"""Запустим процедуру построения"""

# ╔═╡ 0142c5e8-0fdf-11f0-3ac6-4df71525cc7a
md"Визуализируем полученное множество"

# ╔═╡ 4aa35833-7758-49df-a74c-8c14f00958d3
@bind output_filename confirm(@htl("Choose output filename $(TextField())"))

# ╔═╡ 0142c5fc-0fdf-11f0-2342-032b21cc0b09
md"![](results/T2_julia_set_c_0.285_0.000.png)"

# ╔═╡ 0142c624-0fdf-11f0-1886-095a02009662
md"![](results/T2_julia_set_c_neg0.835_neg0.2321.png)"

# ╔═╡ 9365e80f-a23e-4452-b243-56b3feba7621
md"""### Прочее"""

# ╔═╡ a16a215a-9b1d-4d97-8ed8-225780d51013
function iim_configuration_from_ui(method_params)
	return (configuration = InvIterMethodConfiguration(
				method_params.n, #niters 
				full::InversionBranchMethod, #inversion
				 ComplexGridConfiguration(3, 
							-method_params.d - method_params.d*1im,
							 method_params.d + method_params.d*1im), # build_domain 
				),
			init_points  = reduce(vcat, gridmk(ComplexGridConfiguration(1, 
								-method_params.l - method_params.l*1im,
								 method_params.l + method_params.l*im)))
		   )
end

# ╔═╡ 9acfc88c-0a27-4c5d-9e64-1135c0671ffc
begin
	struct SVG
		content
	end

	function Base.show(io::IO, m::MIME"image/svg+xml", s::SVG)
		write(io, s.content)
	end

	SVG
end

# ╔═╡ 0142c644-0fdf-11f0-0c3b-5be86ef27fb0
const im_axes = """<svg xmlns="http://www.w3.org/2000/svg" width="300" height="300" viewBox="0 0 300 300" style="max-width:100%;height:auto;color:black">

<defs>
    <marker id="triangle" viewBox="0 0 10 10" refX="1" refY="5" markerUnits="strokeWidth" markerWidth="10" markerHeight="10" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="currentColor"></path>
    </marker>
  </defs>

<rect width="300" height="300" fill="white" rx="1em"/>

<line x1="10" x2="290" y1="150" y2="150" stroke="currentColor" marker-start="url(#triangle)" marker-end="url(#triangle)"></line>

<line y1="10" y2="290" x1="150" x2="150" stroke="currentColor" marker-start="url(#triangle)" marker-end="url(#triangle)"></line>
<circle cx="150" cy="150" r="3" fill="currentColor"></circle>

<g style="    font-family: math;
    font-style: italic;
    font-size: 17px;">

<g transform="translate(260, 140)">
<text fill="currentColor">Real</text>
</g>

<g transform="translate(170, 50) rotate(-90)">
<text text-anchor="middle" fill="currentColor">Imaginary</text>
</g>


<g transform="translate(150, 150)">
<text text-anchor="end" fill="currentColor" dx="-10" dy="20">0</text>
</g>



<g transform="translate(270, 150)">
<line y1="-5" y2="5" stroke="currentColor"></line>
<text text-anchor="middle" fill="currentColor" dy="20">1</text>
</g>
<g transform="translate(30, 150)">
<line y1="-5" y2="5" stroke="currentColor"></line>
<text text-anchor="middle" fill="currentColor" dy="20">-1</text>
</g>





<g transform="translate(150, 30)">
<line x1="-5" x2="5" stroke="currentColor"></line>
<text text-anchor="end" fill="currentColor" dy=".5ch" dx="-10">i</text>
</g>
<g transform="translate(150, 270)">
<line x1="-5" x2="5" stroke="currentColor"></line>
<text text-anchor="end" fill="currentColor" dy=".5ch" dx="-10">-i</text>
</g>



</g>

</svg>""" |> SVG;

# ╔═╡ 05132c1e-5d8f-4d20-9f12-fd6cd36d15c0
function ComplexNumberPicker(; default::Union{Real,Complex,Nothing}=nothing)
	t(x) = (x - 150.0) / 120.0
	tinv(x) = x * 120.0 + 150.0
	
	default_coord = default === nothing ? nothing : 
		PlutoImageCoordinatePicker.ClickCoordinate(300, 300, 
			tinv(real(default)), 
			tinv(-imag(default))
		)
	
	PlutoUI.Experimental.transformed_value(ImageCoordinatePicker(im_axes; default=default_coord, hint=false)) do point
		if point === nothing
			nothing
		else
			round(t(point.x) - im * t(point.y), digits=3)
		end
	end
end

# ╔═╡ db90c4f7-1b31-4761-922b-f31c38f9f9d0
@bind c confirm(ComplexNumberPicker(default=.285+.0im))

# ╔═╡ 0d825988-e841-40ff-98f6-abaf52c43757
c

# ╔═╡ 0142c58e-0fdf-11f0-3162-a5cfe3d61cba
begin
	z = ComplexVariable()
	system = PolyHolomorphicDS(c + z^2)
end

# ╔═╡ 0142c5de-0fdf-11f0-2b1b-3f8d30832926
js = BuildJuliaSet(system, 
	iim_configuration_from_ui(iim_parameters).configuration, iim_configuration_from_ui(iim_parameters).init_points
);

# ╔═╡ 0142c5f2-0fdf-11f0-315c-e7983fc65004
begin
julia_set = scatter(real(js.elems[js.marks]), imag(js.elems[js.marks]), m=0, ms=0, title="Julia set for f(z) = z² + $(round(c, digits=3))", legend=false);
plot!(size=(1000,1000))
xlabel!("Re")
ylabel!("Im")
end

# ╔═╡ 54687aa4-b3e7-444d-9195-b68768a6a77c
savefig(julia_set, if output_filename != "" "results/"*output_filename else "results/problem-02" end)

# ╔═╡ Cell order:
# ╟─0142c4dc-0fdf-11f0-3403-77788dee9947
# ╟─013db4c4-0fdf-11f0-0aaf-c92fc9d78215
# ╟─dd7de283-d8a8-47ed-96a1-98ee4d715abd
# ╟─0142570c-0fdf-11f0-2ec7-1378fe6665f6
# ╟─0142c5ac-0fdf-11f0-11ec-6919b42c22a9
# ╠═9cd321b6-f105-434c-b0ae-5fc44c6a6988
# ╟─0142c55c-0fdf-11f0-1369-d1e8b81243e6
# ╠═db90c4f7-1b31-4761-922b-f31c38f9f9d0
# ╠═0d825988-e841-40ff-98f6-abaf52c43757
# ╠═0142c58e-0fdf-11f0-3162-a5cfe3d61cba
# ╟─bfc236f3-68ed-4da4-8fcf-8cea6363235d
# ╠═0142c5de-0fdf-11f0-2b1b-3f8d30832926
# ╟─0142c5e8-0fdf-11f0-3ac6-4df71525cc7a
# ╠═0142c5f2-0fdf-11f0-315c-e7983fc65004
# ╟─4aa35833-7758-49df-a74c-8c14f00958d3
# ╟─54687aa4-b3e7-444d-9195-b68768a6a77c
# ╟─0142c5fc-0fdf-11f0-2342-032b21cc0b09
# ╟─0142c624-0fdf-11f0-1886-095a02009662
# ╟─9365e80f-a23e-4452-b243-56b3feba7621
# ╟─a16a215a-9b1d-4d97-8ed8-225780d51013
# ╟─05132c1e-5d8f-4d20-9f12-fd6cd36d15c0
# ╟─9acfc88c-0a27-4c5d-9e64-1135c0671ffc
# ╟─0142c644-0fdf-11f0-0c3b-5be86ef27fb0
