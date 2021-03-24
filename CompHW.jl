### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 8df467aa-8bd9-11eb-088e-218f062da936
md"# Digital Communication Systems
## Computational Homework (TI0069)

*Filipe Pereira de Farias - 397071, UFC*

This notebook presents the implementations of the package DigCommSys, developed for the Digital Communication Systems course at UFC. The codes are in [author's GitHub](https://github.com/filippfarias/DigCommSys) and imported below.
"

# ╔═╡ 25bc3a22-8be4-11eb-1414-0b8e408a7b54
begin
	import Pkg;
	Pkg.add(url="https://github.com/filippfarias/DigCommSys")
end

# ╔═╡ 9544fbf0-8bd9-11eb-2a08-ff3364fb0c58
#using Plots, DigCommSys

# ╔═╡ f1835c24-8bdc-11eb-2e3e-ffaca518a99c
md"""
#### Question 1
###### 1. Constellation Average Energy

To evaluate the average energy for each $M$ constellation of a rectangular QAM [[Proakis]](#ref-proakis) [[Proakis]](#ans-q11) we need to consider first the **outer sum** of two $\sqrt{M}$--PAM signals, being **in-phase** $A_{mi}$ (real axis) and the **quadrature** $A_{mq}$ (imaginary axis) its components. In that way we can define the base band signal as

$s_m(t) = (A_{mi} + jA_{mq})g(t).$

If we expand it using an orthonormal basis such that

$$\phi_1(t)=\sqrt{\frac{2}{\mathcal{E}_g}}g(t)\cos 2\pi f_c t \ \text{ and } \ \phi_2(t)=-\sqrt{\frac{2}{\mathcal{E}_g}}g(t)\sin 2\pi f_c t$$

then we can represent $s_m(t)$ in the bidimensional form, 

$\mathbf{s}_m=\left( A_{mi}\sqrt{\frac{\mathcal{E}_g}{2}}, A_{mq}\sqrt{\frac{\mathcal{E}_g}{2}} \right).$

The energy of $\mathbf{s}_m$ will be given by

$\mathcal{E}_m =||\mathbf{s}_m||^2=\frac{\mathcal{E}_g}{2}(A_{mi}^2+A_{mq}^2).$

In the caso of the rectangular QAM, where $A_m=\pm1,\pm3,\dots,\pm(\sqrt{M}-1)$ for both in-phase and quadrature, we can sum the energy of each symbol at each direction, i.e.

$\begin{aligned}
\mathcal{E}_\text{avg} &= \frac{1}{M} \frac{\mathcal{E}_g}{2}\sum_{m=1}^{\sqrt{M}} \sum_{n=1}^{\sqrt{M}}\left(A_{m}^{2}+A_{n}^{2}\right) \\ 
&= \frac{1}{M} \frac{\mathcal{E}_g}{2}\sum_{m=1}^{\sqrt{M}} \left(\sum_{n=1}^{\sqrt{M}}(A_{m}^{2}+A_{n}^{2})\right)\\
&= \frac{1}{M} \frac{\mathcal{E}_g}{2}\sum_{m=1}^{\sqrt{M}} \left(\sqrt{M}A_{m}^{2}+\sum_{n=1}^{\sqrt{M}}A_{n}^{2}\right) \\
&= \frac{1}{M} \frac{\mathcal{E}_g}{2} \sqrt{M}\left(\sum_{m=1}^{\sqrt{M}}A_{m}^{2}+\sum_{n=1}^{\sqrt{M}}A_{n}^{2}\right)
\end{aligned}$

We can manipulate $A_m=\pm1,\pm3,\dots,\pm(\sqrt{M}-1) = (2m - 1 -\sqrt{M}), \ m=1,\dots,\sqrt{M}$ in the way that the negative amplitudes can be defined as

$A_m^- =-1,-3,\dots,-(\sqrt{M}-1) =-2m-1, \ m = 0,\dots,\sqrt{M}/2-1.$

"""


# ╔═╡ 57549df8-8c3f-11eb-16a5-c33fe4560e41
md"""
Then $A_m = [A_m^-,A_m^+]$ and the positive ones are $A_m^+=-A_m^-=-(-2m-1)$ with $m = 0,\dots,\sqrt{M}/2-1$. Back to the average energy, we have that

$\begin{aligned}
\sum_{m=1}^{\sqrt{M}}A_{m}^{2}+\sum_{n=1}^{\sqrt{M}}A_{n}^{2}&= \left(\sum_{m=0}^{\sqrt{M}/2-1}{A_{m}^+}^{2}+\sum_{m=0}^{\sqrt{M}/2-1}{A_{m}^-}^{2}\right) +
\left(\sum_{n=0}^{\sqrt{M}/2-1}{A_{n}^+}^{2}+\sum_{n=0}^{\sqrt{M}/2-1}{A_{n}^-}^{2}\right) \\
&= \left(\sum_{m=0}^{\sqrt{M}/2-1}{A_{m}^+}^{2}+\sum_{m=0}^{\sqrt{M}/2-1}{A_{m}^-}^{2}\right) +
\left(\sum_{n=0}^{\sqrt{M}/2-1}{A_{n}^+}^{2}+\sum_{n=0}^{\sqrt{M}/2-1}{A_{n}^-}^{2}\right) \\
\end{aligned}$

As $(2m+1)^2 = (-2m-1)^2$ and $m=n, \ m=0,\dots,\sqrt{m}/2-1$, we obtain

$\begin{aligned}
& \left(\sum_{m=0}^{\sqrt{M}/2-1}{A_{m}^+}^{2}+\sum_{m=0}^{\sqrt{M}/2-1}{A_{m}^-}^{2}\right) +
\left(\sum_{n=0}^{\sqrt{M}/2-1}{A_{n}^+}^{2}+\sum_{n=0}^{\sqrt{M}/2-1}{A_{n}^-}^{2}\right) =\\
&= 4\sum_{m=0}^{\sqrt{M}/2-1}{A_{m}^+}^{2} \\
&= 4\sum_{m=0}^{\sqrt{M}/2-1}{(2m+1)}^{2} \\
&= 4\sum_{m=0}^{\sqrt{M}/2-1} (4m^2 + 4m + 1) 
\end{aligned}$

"""

# ╔═╡ 9dbdeb40-8c24-11eb-1569-6f010c7f5ed5
md"""
From another results [[Wiki1]](#ref-wiki1) we know that

$$
 \sum_{m=0}^{\sqrt{M}/2-1}4m^2 = 4\cdot\frac{(\sqrt{M}/2-1)(\sqrt{M}/2)(\sqrt{M}-1)}{6},$$
$$\sum_{m=0}^{\sqrt{M}/2-1}4m = 4\cdot\frac{(\sqrt{M}/2-1)(\sqrt{M}/2)}{2}, \, \sum_{n=0}^{\sqrt{M}/2-1} 1 = \frac{\sqrt{M}}{2}.$$

What leads to

$\begin{aligned}
&4 \left( 4\cdot \frac{(\sqrt{M}/2-1)(\sqrt{M}/2)(\sqrt{M}-1)}{6} + 4\cdot\frac{(\sqrt{M}/2-1)(\sqrt{M}/2)}{2} + \frac{\sqrt{M}}{2} \right) =\\
&= \frac{2\sqrt{M}(M-1)}{3}
\end{aligned}$
"""

# ╔═╡ 98e4d1ba-8c45-11eb-2e90-21ecfb046773
md"""
Finally, we have that 

$\mathcal{E}_\text{avg} =  \frac{1}{M} \frac{\mathcal{E}_g}{2} \sqrt{M}\frac{2\sqrt{M}(M-1)}{3} = \mathcal{E}_g\frac{(M-1)}{3}.$
"""

# ╔═╡ df072ef6-8c63-11eb-0064-8328864e88cf
md"""
###### 2. Minimum distance
"""

# ╔═╡ bd6f75c6-8c63-11eb-1541-f17f0f1803bb
html"""
<h6> 3. Modulator (bit-symbol mapping) </h6><br>

<p>The M-QAM modulator is implemented in <a href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/QAM.jl" target="_blank"><code>QAM.jl</code></a> at the function <code>MQAM</code>. First of all, the function has two arguments, the order of the modulation <code>M</code> and the pulse <code>energy</code>. In the lines<br></p>
"""

# ╔═╡ 45b663b0-8c67-11eb-1302-df7ec6e9452d
md"""
	
	I = ((range(0,2*(m-1),step=d) .- (2*m-2)/2) .+0*1im);
	Q = I.*1im
is evaluated the **in-phase** ```I``` and  **quadrature** amplitudes ```Q``` for the M-QAM constellation, where ```1im``` represents the imaginary $i$. Next we compute the **outer sum** by ```(Q.+transpose(I))``` in line

	constellation = (Q.+transpose(I))/√avgEnergy
obtaining the constellation normalised by the **root of the average energy** ```√avgEnergy```.

	alphabet = GrayCode(Int(log2(M)))
"""

# ╔═╡ 684ef8fc-8c69-11eb-0947-771929a7f355
html"""
<p>The next step is to construct the <code>alphabet</code> according to the <strong>Gray code</strong>. The <code>GrayCode</code> function implemented <a href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/GrayCode.jl" target="_blank">here</a> is based on <a href="#ref-wiki2">[Wiki2]</a>, by recursive concatenating and prefixing. The value inserted in <code>GrayCode</code> is equivalent to the number of bits possible for a M-QAM, i.e. <code>log2(M)</code>. The next steps are needed just to fix the way we construct the arrays and to fill the Gray code correctly. The function returns the <code>constellation</code> and the <code>alphabet</code>.
"""

# ╔═╡ 717edf16-8c6b-11eb-1b71-a51e7d799605
plotly();

# ╔═╡ 0d70be04-8bf6-11eb-0638-a35109c2d36c
md"
---
### References
"

# ╔═╡ 21953290-8bf7-11eb-025f-87f8477b88c2
html"""

<p><a name="ref-proakis">[Proakis]</a> Digital Communication Systems 5th Edition, John G. Proakis and Masoud Salehi.</p>
<p><a name="ref-wiki1">[Wiki1]</a> <a target="_blank" href="https://en.wikipedia.org/wiki/Summation#Powers_and_logarithm_of_arithmetic_progressions">Summation - Powers and logarithm of arithmetic progressions</a>. Wikipedia </p>
<p><a name="ref-wiki2">[Wiki2]</a> <a target="_blank" href="https://en.wikipedia.org/wiki/Gray_code#Constructing_an_n-bit_Gray_code">Gray code - Constructing an n-bit Gray code</a>. Wikipedia </p>
"""

# ╔═╡ 8e37ac06-8c5b-11eb-3ca4-07b42d5b841a
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));

# ╔═╡ bf54e3d6-8c5b-11eb-1cd9-2b7bc38fc404
almost(text) = Markdown.MD(Markdown.Admonition("warning", "Almost there!", [text]));

# ╔═╡ bf554652-8c5b-11eb-30f1-b9352a6f77cf
keep_working(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]));

# ╔═╡ bf5fc848-8c5b-11eb-16e1-b7c9b467c87a
correct(title,text) = Markdown.MD(Markdown.Admonition("correct", title, [text]));

# ╔═╡ 7dfecccc-8c5c-11eb-01e5-1f9b7f06dbc8
begin 
	M_ = [4,16,64];
	e_avg = (M_ .- 1) ./ 3;
	correct("Answer for Question 1.1",md"""
For $M=\{4,16,64\}$, considering $\mathcal{E}_g=1$, the average energy is, respectively, $ $[e_avg[1]], $[e_avg[2]], $[e_avg[3]]. $""");
end

# ╔═╡ 8036d610-8c61-11eb-2cef-43773fdbf97c

begin 
	correct("Answer for Question 1.2",md"""
Being the QAM rectangular, we know that its module is $\sqrt{\mathcal{E}_g}=1$ which is half diagonal of the 4-QAM. Then by Pythagorean theorem we know that $2d^2 = 4\mathcal{E}_g$, what implies that $d_\text{min} = \sqrt{2\mathcal{E}_g}=\sqrt{2}$. By symmetry we can note this minimum distance holds for M-QAM.""");
end

# ╔═╡ 8d99049a-8c6b-11eb-3118-51c16e0927bf
begin 
	correct("Answer for Question 1.3",md"""
Select the modulation: `M = ` $(@bind M html"<select><option value=4>4</option><option value=16>16</option><option value=64>64</option></select>")""");
end

# ╔═╡ 78f40b68-8c6c-11eb-2925-15974d7dccb5
begin
	(alphabetQ13,constellationQ13) = MQAM(parse(Int64,M));
	scatter(real(constellationQ13[:]),imag(constellationQ13[:]),legend=false,lims=(-1.175,1.175));
	annotate!((real(constellationQ13[:]),imag(constellationQ13[:]) .+ .06,alphabetQ13[:],10));
end

# ╔═╡ Cell order:
# ╟─8df467aa-8bd9-11eb-088e-218f062da936
# ╠═25bc3a22-8be4-11eb-1414-0b8e408a7b54
# ╠═9544fbf0-8bd9-11eb-2a08-ff3364fb0c58
# ╟─f1835c24-8bdc-11eb-2e3e-ffaca518a99c
# ╟─57549df8-8c3f-11eb-16a5-c33fe4560e41
# ╟─9dbdeb40-8c24-11eb-1569-6f010c7f5ed5
# ╟─98e4d1ba-8c45-11eb-2e90-21ecfb046773
# ╠═7dfecccc-8c5c-11eb-01e5-1f9b7f06dbc8
# ╟─df072ef6-8c63-11eb-0064-8328864e88cf
# ╟─8036d610-8c61-11eb-2cef-43773fdbf97c
# ╟─bd6f75c6-8c63-11eb-1541-f17f0f1803bb
# ╟─45b663b0-8c67-11eb-1302-df7ec6e9452d
# ╟─684ef8fc-8c69-11eb-0947-771929a7f355
# ╠═717edf16-8c6b-11eb-1b71-a51e7d799605
# ╠═8d99049a-8c6b-11eb-3118-51c16e0927bf
# ╟─78f40b68-8c6c-11eb-2925-15974d7dccb5
# ╟─0d70be04-8bf6-11eb-0638-a35109c2d36c
# ╟─21953290-8bf7-11eb-025f-87f8477b88c2
# ╟─8e37ac06-8c5b-11eb-3ca4-07b42d5b841a
# ╟─bf54e3d6-8c5b-11eb-1cd9-2b7bc38fc404
# ╟─bf554652-8c5b-11eb-30f1-b9352a6f77cf
# ╟─bf5fc848-8c5b-11eb-16e1-b7c9b467c87a
