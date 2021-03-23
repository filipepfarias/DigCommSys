### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 9544fbf0-8bd9-11eb-2a08-ff3364fb0c58
using Plots, DigCommSys

# ╔═╡ 8df467aa-8bd9-11eb-088e-218f062da936
md"# Digital Communication Systems
## Computational Homework (TI0069)

*Filipe Pereira de Farias - 397071, UFC*

This notebook presents the implementations of the package DigCommSys, developed for the Digital Communication Systems course at UFC. The codes are in [author's GitHub](https://github.com/filippfarias/DigCommSys) and imported below.
"

# ╔═╡ 25bc3a22-8be4-11eb-1414-0b8e408a7b54
# begin
# 	import Pkg;
# 	Pkg.add(url="https://github.com/filippfarias/DigCommSys")
# end

# ╔═╡ f1835c24-8bdc-11eb-2e3e-ffaca518a99c
md"""
#### Question 1
1. Constellation Average Energy

To evaluate the average energy for each $M$ constellation of a rectangular QAM [[Proakis]](#ref-proakis) we need to consider first the **outer sum** of two $\sqrt{M}$--PAM signals, being **in-phase** $A_{mi}$ (real axis) and the **quadrature** $A_{mq}$ (imaginary axis) its components. In that way we can define the base band signal as

$s_m(t) = (A_{mi} + jA_{mq})g(t).$

If we expand it using an orthonormal basis such that

$$\phi_1(t)=\sqrt{\frac{2}{\mathcal{E}_g}}g(t)\cos 2\pi f_c t \ \text{ and } \ \phi_2(t)=-\sqrt{\frac{2}{\mathcal{E}_g}}g(t)\sin 2\pi f_c t$$

then we can represent $s_m(t)$ in the bidimensional form, 

$\mathbf{s}_m=\left( A_{mi}\sqrt{\frac{\mathcal{E}_g}{2}}, A_{mq}\sqrt{\frac{\mathcal{E}_g}{2}} \right).$

The energy of $\mathbf{s}_m$ will be given by

$\mathcal{E}_m =||\mathbf{s}_m||^2=\frac{\mathcal{E}_g}{2}(A_{mi}^2+A_{mq}^2).$

In the caso of the rectangular QAM, where $A_m=\pm1,\pm3,\dots,\pm(\sqrt{M}-1)$ for both in-phase and quadrature, we can sum the energy of each symbol at each direction, i.e.

$\begin{aligned}
\mathcal{E}_\text{avg} &= \frac{1}{M} \frac{\mathcal{E}_g}{2}\sum_{m=1}^{\sqrt{M}} \sum_{n=1}^{\sqrt{M}}\left(A_{m}^{2}+A_{n}^{2}\right)
\end{aligned}$

We can manipulate $A_m=\pm1,\pm3,\dots,\pm(\sqrt{M}-1) = (2m - 1 -\sqrt{M}), \ m=1,\dots,\sqrt{M}$ in the way that the negative amplitudes can be defined as

$A_m^- =-1,-3,\dots,-(\sqrt{M}-1) =-2m-1, \ m = 0,\dots,\sqrt{M}/2-1.$

Then the positive ones $A_m^+=-A_m^-=-(-2m-1)$ with $m = 0,\dots,\sqrt{M}/2-1$. Back to the average energy, we have that

$\begin{aligned}
\sum_{m=1}^{\sqrt{M}} \sum_{n=1}^{\sqrt{M}}\left(A_{m}^{2}+A_{n}^{2}\right) &= \sum_{m=0}^{\sqrt{M}/2-1} \sum_{n=0}^{\sqrt{M}/2-1}\left({A_{m}^-}^{2}+{A_{n}^-}^{2}\right) + \sum_{m=0}^{\sqrt{M}/2-1} \sum_{n=0}^{\sqrt{M}/2-1}\left((-A_{m}^+)^{2}+(-A_{n}^+)^{2}\right) \\
&=\sum_{m=0}^{\sqrt{M}/2-1} \sum_{n=0}^{\sqrt{M}/2-1}\left({(-2m-1)}^{2}+{(-2n-1)}^{2}\right) + \\
& \quad + \sum_{m=0}^{\sqrt{M}/2-1} \sum_{n=0}^{\sqrt{M}/2-1}\left((2m+1)^{2}+(2n+1)^{2}\right) \\
& = 2\sum_{m=0}^{\sqrt{M}/2-1} \sum_{n=0}^{\sqrt{M}/2-1}\left((2m+1)^{2}+(2n+1)^{2}\right) \\ 
& = 2\sum_{m=0}^{\sqrt{M}/2-1} \sum_{n=0}^{\sqrt{M}/2-1}\left(4m^2 + 4n^2 + 4m + 4n + 2\right) \\ 
& = 2\sum_{m=0}^{\sqrt{M}/2-1}\left( \sum_{n=0}^{\sqrt{M}/2-1}4m^2 + \sum_{n=0}^{\sqrt{M}/2-1}4n^2\right. \\
& \quad \left.+ \sum_{n=0}^{\sqrt{M}/2-1}4m + \sum_{n=0}^{\sqrt{M}/2-1}4n + \sum_{n=0}^{\sqrt{M}/2-1}2\right) \\ 
\end{aligned}$

"""


# ╔═╡ 9dbdeb40-8c24-11eb-1569-6f010c7f5ed5
md"""
From another results [[Wiki1]](#ref-wiki1) we know that

$\begin{aligned}
 \sum_{n=0}^{\sqrt{M}/2-1}4m^2 = 4m^2 \cdot \sqrt{M}/2, \quad \sum_{n=0}^{\sqrt{M}/2-1}4n^2 = 4\cdot\frac{(\sqrt{M}/2-1)(\sqrt{M}/2)(\sqrt{M}+1)}{6} \\ \sum_{n=0}^{\sqrt{M}/2-1}4m = 4m \cdot \sqrt{M}/2, \, \sum_{n=0}^{\sqrt{M}/2-1}4n = 4\cdot\frac{(\sqrt{M}/2-1)(\sqrt{M}/2)}{6}, \, \sum_{n=0}^{\sqrt{M}/2-1}2 = 2 \cdot  \sqrt{M}/2
\end{aligned}$
"""

# ╔═╡ df77ee44-8c29-11eb-2b3d-ffa5dfc6a409
md"""
Then we substitute

$\begin{aligned}
& 2 \sum_{m=0}^{\sqrt{M}/2-1}\left( 4m^2 \cdot \sqrt{M}/2 + 4\cdot\frac{(\sqrt{M}/2-1)(\sqrt{M}/2)(\sqrt{M}+1)}{6}\right. \\
& \quad \left.+4m \cdot \sqrt{M}/2 + 4\cdot\frac{(\sqrt{M}/2-1)(\sqrt{M}/2)}{2} +  2 \cdot  \sqrt{M}/2\right) \\
&= 2 \cdot 4 \frac{(\sqrt{M}/2-1)(\sqrt{M}/2)(\sqrt{M}+1)}{6}\cdot \sqrt{M}/2 + \\
& \quad + 2\cdot 4\cdot\frac{(\sqrt{M}/2-1)(\sqrt{M}/2)(\sqrt{M}+1)}{6} \cdot \sqrt{M}/2+\\
& \quad + 2 \cdot 4\cdot\frac{(\sqrt{M}/2-1)(\sqrt{M}/2)}{2} \cdot \sqrt{M}/2 + \\
& \quad + 2 \cdot 4\cdot\frac{(\sqrt{M}/2-1)(\sqrt{M}/2)}{2} \cdot \sqrt{M}/2 + \\
& \quad + 2\cdot 2 \cdot \sqrt{M}/2 \cdot \sqrt{M}/2 \\
&= 2 M \frac{(\sqrt{M}/2-1)(\sqrt{M}+1)}{6} + 2M\cdot\frac{(\sqrt{M}/2-1)(\sqrt{M}+1)}{6}+\\
& \quad + 2 M\cdot\frac{(\sqrt{M}/2-1)}{2} +  2 M\cdot\frac{(\sqrt{M}/2-1)}{2} + M \\
&=
\end{aligned}$
"""

# ╔═╡ 6f445cfe-8be6-11eb-14da-f11e24d94094
(alphabet,constellation) = MQAM(64)

# ╔═╡ c5273a9a-8bd9-11eb-134f-f34babe7faa2
begin
	scatter(real(constellation[:]),imag(constellation[:]));
	annotate!((real(constellation[:]),imag(constellation[:]).+.1,alphabet[:]))
end

# ╔═╡ 0e7f5c2a-8bff-11eb-1b1e-839b482632e5
sum(abs.(constellation[:]))/64

# ╔═╡ 0d70be04-8bf6-11eb-0638-a35109c2d36c
md"
---
### References
"

# ╔═╡ 21953290-8bf7-11eb-025f-87f8477b88c2
html"""

<span name="ref-proakis">[Proakis] Digital Communication Systems 5th Edition, John G. Proakis and Masoud Salehi.</span><br>
<span name="ref-wiki1">[Wiki1] <a target="_blank" href="https://en.wikipedia.org/wiki/Summation#Powers_and_logarithm_of_arithmetic_progressions">Summation - Powers and logarithm of arithmetic progressions</a>. Wikipedia </span>
"""

# ╔═╡ Cell order:
# ╠═8df467aa-8bd9-11eb-088e-218f062da936
# ╠═25bc3a22-8be4-11eb-1414-0b8e408a7b54
# ╠═9544fbf0-8bd9-11eb-2a08-ff3364fb0c58
# ╟─f1835c24-8bdc-11eb-2e3e-ffaca518a99c
# ╠═9dbdeb40-8c24-11eb-1569-6f010c7f5ed5
# ╠═df77ee44-8c29-11eb-2b3d-ffa5dfc6a409
# ╠═6f445cfe-8be6-11eb-14da-f11e24d94094
# ╠═c5273a9a-8bd9-11eb-134f-f34babe7faa2
# ╠═0e7f5c2a-8bff-11eb-1b1e-839b482632e5
# ╟─0d70be04-8bf6-11eb-0638-a35109c2d36c
# ╟─21953290-8bf7-11eb-025f-87f8477b88c2
