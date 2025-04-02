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

# ╔═╡ 040bc779-2a7d-43af-82d2-2608f02bdc60
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
	using Polynomials
end

# ╔═╡ 96058834-100a-11f0-0a28-379ac44bb829
md"# Задача III. Исследование динамической системы"

# ╔═╡ 960a0710-100a-11f0-2f62-ef6acc36b76d
md"""
Рассматривается динамическая система порождаемая отображением
$$z \mapsto R(z) = \lambda z + z^2$$
где $\lambda = \exp(2\pi i \alpha), \alpha = \frac{1}{20}$. 


Необходимо найти неподвижные точки, определить тип их устойчивости, цикл периода 2 и построить приближение к инвариантному множеству
"""

# ╔═╡ 960a6ffc-100a-11f0-228d-156a0cd0a72a
begin
	α = 1/20
	λ = exp(2*1im*π*α)
	
	z = ComplexVariable()
	R = λ * z + z^2
end

# ╔═╡ 960a702e-100a-11f0-265d-2b8b91b8b209
md"### Неподвижные точки и их тип устойчивости"

# ╔═╡ 960a7042-100a-11f0-2cd6-1b5df6b2bce1
md"Неподвижные точки системы удовлетворяют уравнению $R(z) = z$:"

# ╔═╡ 960a704c-100a-11f0-2a55-03a9d1844b3d
stability_points = roots(R - z)

# ╔═╡ 960a7058-100a-11f0-3e36-3584c5a2e52e
md"Определим тип устойчивости этих точек. Точка $z_0 = 0$ будет нейтральной параболической неподвижной точкой, так как производная имеет вид $R'(0) = \lambda = \exp(2 \pi i \alpha)$, где $\alpha$ рациональное. Вторая неподвижная точка  $z_1 \approx 0.0489 -  0.3090 \cdot i$ будет отталкивающей неподвижной точкой так как $|R'(z_1)| > 1$ (см. вычисление модуля производной ниже)"

# ╔═╡ 960a706a-100a-11f0-1b1d-4bfa020097ee
"|R'(z)| for stability points: " * string(abs.(derivative(R)(stability_points)))

# ╔═╡ 960a707e-100a-11f0-0e01-25654efb24a0
md"### Циклы периода 2"

# ╔═╡ 960a7092-100a-11f0-22e0-8d1ca1f35ac3
md"""
После двух итераций отображения, точка $z$ переходит в 

$$z \mapsto R(R(z)) = (\lambda z + z^2) \cdot (\lambda + \lambda z + z^2) = R(z) \cdot (\lambda + R(z))$$
"""

# ╔═╡ 960a709c-100a-11f0-359d-4715e5c979a6
second_iterration_R = R * (λ + R);

# ╔═╡ 960a70a6-100a-11f0-04c5-99993ac471f1
md"Чтобы найти циклы периода 2 надо численно решить уравнение $R(R(z)) = z$ и взять те корни, для которых $R(z) \neq z$:"

# ╔═╡ 960a70b0-100a-11f0-2c7a-c7db00f17212
begin
	cycles = roots(second_iterration_R - z)
	cycles = filter(z -> minimum(abs.(stability_points .- z)) > 0.00001, cycles)
	cycles
end

# ╔═╡ 960a70b8-100a-11f0-0e34-d9ace498ef97
md"### Построение инвариантного множества"

# ╔═╡ 960a70c4-100a-11f0-2d01-09637414c60d
md"Используем метод обратных итераций и построим множество Жюлиа для данной системы. (описание параметров дано в тетради `problem-02.jl`)"

# ╔═╡ db1a12d4-c7dc-4868-bf67-b75d3a017309
system = PolyHolomorphicDS(R);

# ╔═╡ a91c44b5-dac3-4435-9e14-c17ea2c49744
@bind iim_parameters confirm(
	combine() do Child
	@htl("""
		<h5>Параметры метода</h5>
		<ul>
		$([
			@htl("<li> Параметр границы области, d = $(Child("d", NumberField(0:0.1:3; default=2.0)))"),
			@htl("<li>Параметры начальной сетки точек. Параметр границы сетки, l = $(Child("l", NumberField(0:0.1:2; default=1))), шаг h =  $(Child("h", NumberField(0:0.05:1; default=0.1)))"),
			@htl("<li>Число итераций: $(Child("n", NumberField(0:50; default=40)))")
		])
		</ul>
		""")
	end	
	)

# ╔═╡ 745a7251-0410-4df7-a729-9eddf531a984
@bind output_filename confirm(@htl("Choose output filename $(TextField())"))

# ╔═╡ 960a7100-100a-11f0-1774-89e627767414
md"![](results/T3_julia_set.png)"

# ╔═╡ 19d52f5b-6af7-47bf-95b4-33027dd8582b
md"""# Прочее"""

# ╔═╡ c0c26cdc-b304-4aa1-8965-987ba0a85880
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

# ╔═╡ 960a70ec-100a-11f0-1bc3-9f67380920c2
js = BuildJuliaSet(system, 
	iim_configuration_from_ui(iim_parameters).configuration, iim_configuration_from_ui(iim_parameters).init_points
);

# ╔═╡ 960a70f6-100a-11f0-062e-49322a0cb282
begin
	julia_set = scatter(real(js.elems[js.marks]), imag(js.elems[js.marks]), m=0, ms=0, title="Julia set for f(z) = z² + λz", legend=false);
	xlabel!("Re")
	ylabel!("Im")
end

# ╔═╡ 82609208-65b4-405d-869d-de30f50fd57d
savefig(julia_set, "results/"*output_filename)

# ╔═╡ Cell order:
# ╟─040bc779-2a7d-43af-82d2-2608f02bdc60
# ╟─96058834-100a-11f0-0a28-379ac44bb829
# ╟─960a0710-100a-11f0-2f62-ef6acc36b76d
# ╠═960a6ffc-100a-11f0-228d-156a0cd0a72a
# ╟─960a702e-100a-11f0-265d-2b8b91b8b209
# ╟─960a7042-100a-11f0-2cd6-1b5df6b2bce1
# ╟─960a704c-100a-11f0-2a55-03a9d1844b3d
# ╟─960a7058-100a-11f0-3e36-3584c5a2e52e
# ╟─960a706a-100a-11f0-1b1d-4bfa020097ee
# ╟─960a707e-100a-11f0-0e01-25654efb24a0
# ╟─960a7092-100a-11f0-22e0-8d1ca1f35ac3
# ╠═960a709c-100a-11f0-359d-4715e5c979a6
# ╟─960a70a6-100a-11f0-04c5-99993ac471f1
# ╠═960a70b0-100a-11f0-2c7a-c7db00f17212
# ╟─960a70b8-100a-11f0-0e34-d9ace498ef97
# ╟─960a70c4-100a-11f0-2d01-09637414c60d
# ╠═db1a12d4-c7dc-4868-bf67-b75d3a017309
# ╟─a91c44b5-dac3-4435-9e14-c17ea2c49744
# ╠═960a70ec-100a-11f0-1bc3-9f67380920c2
# ╟─960a70f6-100a-11f0-062e-49322a0cb282
# ╟─745a7251-0410-4df7-a729-9eddf531a984
# ╟─82609208-65b4-405d-869d-de30f50fd57d
# ╟─960a7100-100a-11f0-1774-89e627767414
# ╟─19d52f5b-6af7-47bf-95b4-33027dd8582b
# ╟─c0c26cdc-b304-4aa1-8965-987ba0a85880
