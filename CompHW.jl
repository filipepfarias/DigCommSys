### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 9544fbf0-8bd9-11eb-2a08-ff3364fb0c58
using Plots, DigCommSys

# ╔═╡ 8df467aa-8bd9-11eb-088e-218f062da936
md"# Digital Communication Systems
## Computational Homework (TI0069)

This notebook presents the implementations of the package DigCommSys, developed for the Digital Communication Systems course at UFC. The codes are in [author's GitHub](https://github.com/filippfarias/DigCommSys) and imported below.
"

# ╔═╡ 25bc3a22-8be4-11eb-1414-0b8e408a7b54
# begin
# 	import Pkg;
# 	Pkg.add(url="https://github.com/filippfarias/DigCommSys")
# end

# ╔═╡ f1835c24-8bdc-11eb-2e3e-ffaca518a99c
md"
#### Question 1
1. Constellation Average Energy

To evaluate the average energy for each $M$ constellation of a rectangular QAM [[Proakis]](#ref-proakis) we need to consider first 

$$A_m = \pm1,\pm3,\dots,\pm(\sqrt{M}-1).$$

This represents the amplitudes of the digital PAM. With a **outer sum** of the arrays containing the **in-phase** $A_{mi}=[\pm1,\pm3,\dots,\pm(\sqrt{M}-1)]$ and the **quadrature** $A_{mq}=[\pm1j,\pm3j,\dots,$ $\pm(\sqrt{M}-1)j]$ components.

"


# ╔═╡ 6f445cfe-8be6-11eb-14da-f11e24d94094
(alphabet,constellation) = MQAM(4)

# ╔═╡ c5273a9a-8bd9-11eb-134f-f34babe7faa2
begin
	scatter(real(constellation[:]),imag(constellation[:]));
	annotate!((real(constellation[:]),imag(constellation[:]).+.1,alphabet[:]))
end

# ╔═╡ 0d70be04-8bf6-11eb-0638-a35109c2d36c
md"
---
### References
"

# ╔═╡ 21953290-8bf7-11eb-025f-87f8477b88c2
html"""

<a name="ref-proakis">[Proakis] Digital Communication Systems 5th Edition, John G. Proakis and Masoud Salehi.</a>

"""

# ╔═╡ Cell order:
# ╟─8df467aa-8bd9-11eb-088e-218f062da936
# ╠═25bc3a22-8be4-11eb-1414-0b8e408a7b54
# ╠═9544fbf0-8bd9-11eb-2a08-ff3364fb0c58
# ╠═f1835c24-8bdc-11eb-2e3e-ffaca518a99c
# ╠═6f445cfe-8be6-11eb-14da-f11e24d94094
# ╠═c5273a9a-8bd9-11eb-134f-f34babe7faa2
# ╟─0d70be04-8bf6-11eb-0638-a35109c2d36c
# ╟─21953290-8bf7-11eb-025f-87f8477b88c2
