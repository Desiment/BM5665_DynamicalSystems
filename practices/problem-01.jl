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

# ╔═╡ 6d74a1d4-1af5-47f3-8d55-e57244fe50c4
# ╠═╡ show_logs = false
begin
	import Pkg; Pkg.activate("..") #use project's manifest
	using PlutoUI, HypertextLiteral, Plots; #imports
	import PlutoUI: combine #for method configuration button
	include("../src/SegmentIterationMethod.jl") # method implementation
end

# ╔═╡ 2ea660b7-a2f0-46d5-974e-4866cdbffb3e
md"""
# Практики. Задача 1
Построить приближения к инвариантному множеству дискретной системы 2 порядка методом итераций отрезка (отображения можно брать из списка программы WInSet)
"""

# ╔═╡ d03eed77-0cac-411c-854f-2ef79f2d78f4
md"""
## Решение

Для данной задачи были выбраны следующие системы:
- Система Каталы
- Система Эно

Система Каталы задается следующим отображением:

$$\begin{align*}
	& x_{n+1} = p_1 \cdot x_{n} + y_{n} \\
	&y_{n+1} = p_2 + x_{n}^2
\end{align*}$$

Система Эно задается отображением:

$$\begin{align*}
	& x_{n+1} = 1 + y_{n}  + s_{1}x_n^2 \\
	& y_{n+1} = s_2 x_{n}
\end{align*}$$

"""

# ╔═╡ 188550a0-870c-4f35-bfb4-0d1e5ea15d5b
function cathala(u::Point, parameters::Vector{Float64})
    return (parameters[1] * u[1] + u[2], parameters[2] + u[1] ^ 2)
end

# ╔═╡ 8000d58a-fd65-4ca9-a3c7-fe054c38a365
function henon(u::Point, parameters::Vector{Float64})
    return (1 - parameters[1] * u[1]^2 + u[2], parameters[2] * u[1])
end

# ╔═╡ 4c6081cf-eefe-46e1-8e92-013c1e847f5c
md"""
### Описание процедуры

Рассматриваем динамическую систему, порождаемую отображением $f \colon \mathbb R^2 \to \mathbb R^2$.

Ломаная $L$ состоящая из звеньев $[z_1, z_2], \ldots, [z_{k-1}, z_k], z_i \in \mathbb R^2$ обозначается как $(z_1, \ldots, z_k)$. Под ломанной $f(L)$ понимается ломанная $(f(z_1), \ldots, f(z_k))$. Под суммой ломанных $L_1$ и $L_2$,  $L_1 = (z_1, \ldots z_k)$, $L_2 = (w_1, \ldots, w_m)$ понимается ломанная 

$$L_1 \oplus L_2 = \begin{cases}(z_1, \ldots, z_k, w_1, \ldots w_m) & z_k \neq w_1 \\ (z_1, \ldots, z_k, w_2, \ldots w_m)  & z_k = w_1 \end{cases}$$

Алгоритм построения приближения к инвариантному многообразию в области $D$:
1. Начинаем с некоторой кусочной-ломанной $L_0 = (z_1, \ldots, z_k) = (z_1, z_2) \oplus (z_2, z_3) \oplus \ldots \oplus (z_{k-1}, z_{k})$
2. Кривая $\hat{L}_{i}$ определяется как $\mathrm{Fragmentize}(L_i)$, где $\mathrm{Fragmentize}(L)$ это рекурсивная процедура:
    * Если $L = (z_1, z_2)$ (т.е. $L$ это ломанная из одного звена), причем  $f(L) \subset D$ и $\mathrm{dist}(f(z_1), f(z_2)) \leqslant h$, то $\mathrm{Fragmentize}(L) = (z_1, z_2)$
    * Если $L = (z_1, z_2)$ (т.е. $L$ это ломанная из одного звена), причем $f(L) \subset D$ и $\mathrm{dist}(f(z_1), f(z_2)) > h$, то $L$ делится на две равные части $L = L_{1} \oplus L_2$, где $L_1 = (z_1, z_{\mathrm{mid}}), L_2 = (z_{\mathrm{mid}}, z_2)$ и точка $z_{\mathrm{mid}}$ это середина отрезка $(z_1, z_2)$. В этом случае $\mathrm{Fragmentize}(L) =  \mathrm{Fragmentize}(L_1) \oplus \mathrm{Fragmentize}(L_2)$ 
    * Если $L = (z_1, z_2)$, (т.е. $L$ это ломанная из одного звена) и $f(L) \not\subset D$ то cегмент исключается из рассмотрения: $\mathrm{Fragmentize}(L) = \varnothing$
    * Если $L = L_1 \oplus \ldots \oplus L_k$, то процедура применяется к каждому звену ломанной $\mathrm{Fragmentize}(L) =  \mathrm{Fragmentize}(L_1)  \oplus \ldots \oplus  \mathrm{Fragmentize}(L_k)$
3. Полагаем $L_{i+1} = f(\hat{L}_i), i = 0, \ldots, n$

Дополнительные предположения
- В реализации область $D$ это прямоугольник с центром в $d_0$, шириной $W$ и высотой $H$.
- Стартовая ломанная $L_0$ это всегда некоторый квадрат с центром в $d_0$ и стороной $S$.

Таким образом параметры метода это $h \in \mathbb R_{+}, d_0 \in \mathbb{R}^2, H \in \mathbb R_{+} , W  \in \mathbb R_{+}, S  \in \mathbb R_{+}, n_{\mathrm{iterr}} \in \mathbb{N}$
"""

# ╔═╡ 92f0e28d-95a8-4c11-93bb-78cf12446e94
md"""
### Система Каталы
"""

# ╔═╡ 7423c5e5-55af-4972-a1bd-74e054dfc5b4
@bind cathala_parameters confirm(
	combine() do Child
	@htl("""
		<h5>Параметры системы</h5>
		<ul>
		$([
			@htl("<li>p₁: $(Child("p_1", NumberField(-1:0.01:1; default=0.7)))"),
			@htl("<li>p₂: $(Child("p_2", NumberField(-1:0.01:1; default=-0.82)))")
		])
		</ul>
		""")
	end	
	)

# ╔═╡ e4486237-377d-4e5f-a65e-71f174a6c7dd
@bind cathala_method_parameters confirm(
	combine() do Child
	@htl("""
		<h5>Параметры метода</h5>
		<ul>
		$([
			@htl("<li>Центральная точка d₀. d₀.x = $(Child("d0_x", NumberField(-5:0.1:5; default=0))); d₀.y = $(Child("d0_y", NumberField(-5:0.1:5; default=0)))"),
			@htl("<li>Параметры области D. Высота = $(Child("H", NumberField(0:0.1:10; default=4))); Ширина = $(Child("W", NumberField(0:0.1:10; default=4)))"),
			@htl("<li>Сторона начального квадрата L₀: $(Child("S", NumberField(0:0.1:2; default=1.5)))"),
			@htl("<li>h: $(Child("h", NumberField(0:0.0001:0.1; default=0.001)))"),
			@htl("<li>Число итераций: $(Child("n", NumberField(0:50; default=50)))")
		])
		</ul>
		""")
	end	
	)

# ╔═╡ 4285861d-1132-4049-b076-2787be315925
@bind output_filename_cathala confirm(@htl("Choose output filename $(TextField())"))

# ╔═╡ c5886b94-cc49-4537-8257-369ab176de13
md"""
#### Отображение Эно
"""

# ╔═╡ 67362432-09ac-44e8-b2ce-92e42cfaba70
@bind henon_parameters confirm(
	combine() do Child
	@htl("""
		<h5>Параметры системы</h5>
		<ul>
		$([
			@htl("<li>s₁: $(Child("s_1", NumberField(-1:0.01:1; default=1.4)))"),
			@htl("<li>s₂: $(Child("s_2", NumberField(-1:0.01:1; default=0.3)))")
		])
		</ul>
		""")
	end	
	)

# ╔═╡ adda50bc-5fef-434d-bad8-b1863ca6f341
@bind henon_method_parameters confirm(
	combine() do Child
	@htl("""
		<h5>Параметры метода</h5>
		<ul>
		$([
			@htl("<li>Центральная точка d₀. d₀.x = $(Child("d0_x", NumberField(-5:0.1:5; default=0))); d₀.y = $(Child("d0_y", NumberField(-5:0.1:5; default=0)))"),
			@htl("<li>Параметры области D. Высота = $(Child("H", NumberField(0:0.1:10; default=4))); Ширина = $(Child("W", NumberField(0:0.1:10; default=4)))"),
			@htl("<li>Сторона начального квадрата L₀: $(Child("S", NumberField(0:0.1:2; default=0.75)))"),
			@htl("<li>h: $(Child("h", NumberField(0:0.0001:0.1; default=0.01)))"),
			@htl("<li>Число итераций: $(Child("n", NumberField(0:50; default=20)))")
		])
		</ul>
		""")
	end	
	)

# ╔═╡ d570c279-328e-4e84-a68c-6c949f08e141
@bind output_filename_henon confirm(@htl("Choose output filename $(TextField())"))

# ╔═╡ 9d372eca-6e28-482c-9566-0ae31187d7d3
md"""
# Прочее
"""

# ╔═╡ 2505b7a7-bc17-4504-befc-69a4a26cce74
function segment_iteration_method_configuration_from_ui(method_params)
	d0 = Point((method_params.d0_x,  method_params.d0_y))
	dD = Point((method_params.W / 2, method_params.H / 2))
	return SegmentIterrationConfiguration(
	Pair(d0 .- dD, d0 .+ dD), #domain
	method_params.n, #nitters
	method_params.h, #fragm_dist
	[d0 .+ Point(((sign[1] * method_params.S / 2, sign[2] *  method_params.S / 2)))  for sign in [(-1, -1), (-1, 1), (1, 1), (1, -1), (-1, -1)]] #chain_base
	);
end


# ╔═╡ 3816d364-bfde-46b5-ba3c-a8ad08597a36
s_cathala = SegmentIterration(
	x -> cathala(x,[cathala_parameters.p_1, cathala_parameters.p_2]), segment_iteration_method_configuration_from_ui(cathala_method_parameters)
);

# ╔═╡ 41ea144f-9a7d-4c27-881f-41acec4f87ec
cathala_invariant_set = plot(map(x->x[1], s_cathala), map(x->x[2], s_cathala), label = nothing)

# ╔═╡ 8fd68b8c-ed26-410a-84b0-353c817a74ad
savefig(cathala_invariant_set,  if output_filename_cathala != "" "results/"*output_filename_cathala else "results/problem-01-cathala.png" end);

# ╔═╡ 3cf8e54b-44bd-4e03-8993-ea7f06c0bb06
s_henon = SegmentIterration(
	x -> henon(x,[henon_parameters.s_1, henon_parameters.s_2]), segment_iteration_method_configuration_from_ui(henon_method_parameters)
);

# ╔═╡ af42638e-702d-4a75-868d-e379ce4873e4
henon_invariant_set = plot(map(x->x[1], s_henon), map(x->x[2], s_henon), label = nothing)

# ╔═╡ 70dbeb84-631e-4186-85c2-7037cc3ee1c1
savefig(henon_invariant_set, if output_filename_henon != "" "results/"*output_filename_henon else "results/problem-01-henon.png" end);

# ╔═╡ Cell order:
# ╠═6d74a1d4-1af5-47f3-8d55-e57244fe50c4
# ╟─2ea660b7-a2f0-46d5-974e-4866cdbffb3e
# ╟─d03eed77-0cac-411c-854f-2ef79f2d78f4
# ╠═188550a0-870c-4f35-bfb4-0d1e5ea15d5b
# ╠═8000d58a-fd65-4ca9-a3c7-fe054c38a365
# ╟─4c6081cf-eefe-46e1-8e92-013c1e847f5c
# ╟─92f0e28d-95a8-4c11-93bb-78cf12446e94
# ╟─7423c5e5-55af-4972-a1bd-74e054dfc5b4
# ╟─e4486237-377d-4e5f-a65e-71f174a6c7dd
# ╠═3816d364-bfde-46b5-ba3c-a8ad08597a36
# ╠═41ea144f-9a7d-4c27-881f-41acec4f87ec
# ╟─4285861d-1132-4049-b076-2787be315925
# ╟─8fd68b8c-ed26-410a-84b0-353c817a74ad
# ╟─c5886b94-cc49-4537-8257-369ab176de13
# ╟─67362432-09ac-44e8-b2ce-92e42cfaba70
# ╠═adda50bc-5fef-434d-bad8-b1863ca6f341
# ╠═3cf8e54b-44bd-4e03-8993-ea7f06c0bb06
# ╠═af42638e-702d-4a75-868d-e379ce4873e4
# ╟─d570c279-328e-4e84-a68c-6c949f08e141
# ╟─70dbeb84-631e-4186-85c2-7037cc3ee1c1
# ╟─9d372eca-6e28-482c-9566-0ae31187d7d3
# ╟─2505b7a7-bc17-4504-befc-69a4a26cce74
