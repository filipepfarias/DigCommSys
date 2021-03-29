### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 2780cb8e-901c-11eb-3bf6-6149f19d97c3
using Plots, DigCommSys

# ╔═╡ 8df467aa-8bd9-11eb-088e-218f062da936
md"# Digital Communication Systems
## Computational Homework (TI0069)

*Filipe Pereira de Farias - 397071, UFC*

This notebook presents the implementations of the package DigCommSys, developed for the Digital Communication Systems course at UFC. The codes are in [author's GitHub](https://github.com/filippfarias/DigCommSys) and imported below.
"

# ╔═╡ 25bc3a22-8be4-11eb-1414-0b8e408a7b54
begin
	import Pkg
	Pkg.add(url="https://github.com/filippfarias/DigCommSys")
end

# ╔═╡ f1835c24-8bdc-11eb-2e3e-ffaca518a99c
md"""
#### Question 1
###### 1. Constellation Average Energy

To evaluate the average energy for each $M$ constellation of a rectangular QAM [[Proakis]](#ref-proakis) we need to consider first the **outer sum** of two $\sqrt{M}$--PAM signals, being **in-phase** $A_{mi}$ (real axis) and the **quadrature** $A_{mq}$ (imaginary axis) its components. In that way we can define the base band signal as

$s_m(t) = (A_{mi} + jA_{mq})g(t).$

If we expand it using an orthonormal basis such that

$$\phi_1(t)=\sqrt{\frac{2}{\mathcal{E}_g}}g(t)\cos 2\pi f_c t \ \text{ and } \ \phi_2(t)=-\sqrt{\frac{2}{\mathcal{E}_g}}g(t)\sin 2\pi f_c t$$

then we can represent $s_m(t)$ in the bidimensional form, 

$\mathbf{s}_m=\left( A_{mi}\sqrt{\frac{\mathcal{E}_g}{2}}, A_{mq}\sqrt{\frac{\mathcal{E}_g}{2}} \right).$

The energy of $\mathbf{s}_m$ will be given by

$\mathcal{E}_m =||\mathbf{s}_m||^2=\frac{\mathcal{E}_g}{2}(A_{mi}^2+A_{mq}^2).$

In the case of the rectangular QAM, where $A_m=\pm1,\pm3,\dots,\pm(\sqrt{M}-1)$ for both in-phase and quadrature, we can sum the energy of each symbol at each direction, i.e.

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

# ╔═╡ 53d7bd68-8f3a-11eb-251c-d9b4a34b8a44
md"""
Note that, for **normalized signals**, eu must divide the sent symbols by the **constellation average energy** in order to keep the average energy of the sent signal unitary.
"""

# ╔═╡ bd6f75c6-8c63-11eb-1541-f17f0f1803bb
html"""
<h6> 3. Modulator (bit-symbol mapping) </h6><br>

<p>The M-QAM modulator is implemented in <a href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/QAM.jl" target="_blank"><code>QAM.jl</code></a> at the function <code>MQAM</code>. First of all, the function has two arguments, the order of the modulation <code>M</code> and the pulse <code>energy</code>. In the lines<br></p>
"""

# ╔═╡ 45b663b0-8c67-11eb-1302-df7ec6e9452d
md"""
	
	I = ((2*(1:m) .- 1 .- m)*d/2) .+0*1im;
	Q = I.*1im
is evaluated the **in-phase** ```I``` and  **quadrature** amplitudes ```Q``` for the M-QAM constellation, where ```1im``` represents the imaginary $i$, according to $(2m - 1 -\sqrt{M})$. Next we compute the **outer sum** by ```(Q.+transpose(I))``` in line

	constellation = (Q.+transpose(I))
"""

# ╔═╡ 684ef8fc-8c69-11eb-0947-771929a7f355
html"""
<p>The next step is to construct the <code>alphabet</code> according to the <strong>Gray code</strong>. The <code>GrayCode</code> function implemented <a href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/GrayCode.jl" target="_blank">here</a> is based on <a href="#ref-wiki2">[Wiki2]</a>, by recursive concatenating and prefixing. The value inserted in <code>GrayCode</code> is equivalent to the number of bits possible for a M-QAM, i.e. <code>log2(M)</code>. The next steps are needed just to fix the way we construct the arrays and to fill the Gray code correctly. The function returns the <code>constellation</code> and the <code>alphabet</code>.
"""

# ╔═╡ 4d3d96b6-8ce0-11eb-21c2-ebfbda60351f
md"""
	alphabet = GrayCode(Int(log2(M)))
"""

# ╔═╡ 78f40b68-8c6c-11eb-2925-15974d7dccb5
begin
	plotly();
	_,alphabetQ13,constellationQ13 = MQAM(4,1,unitAveragePower=false);
	scatter(real(constellationQ13[:]),imag(constellationQ13[:]),legend=false,extra_plot_kwargs = KW(:yaxis => KW(:autorange => true),:xaxis => KW(:autorange => true)));
	annotate!((real(constellationQ13[:]),imag(constellationQ13[:]) .+ .06,alphabetQ13[:],10),title="4-QAM");
end

# ╔═╡ 671b4d78-9035-11eb-0a9e-a3abea476a77
begin
	plotly();
	_1,alphabetQ131,constellationQ131 = MQAM(16,1,unitAveragePower=false);
	scatter(real(constellationQ131[:]),imag(constellationQ131[:]),legend=false,extra_plot_kwargs = KW(:yaxis => KW(:autorange => true),:xaxis => KW(:autorange => true)));
	annotate!((real(constellationQ131[:]),imag(constellationQ131[:]) .+ .06,alphabetQ131[:],10),title="16-QAM");
end

# ╔═╡ 68a8d124-9035-11eb-0d05-97e7e4454d10
begin
	plotly();
	_2,alphabetQ132,constellationQ132 = MQAM(64,1,unitAveragePower=false);
	scatter(real(constellationQ132[:]),imag(constellationQ132[:]),legend=false,extra_plot_kwargs = KW(:yaxis => KW(:autorange => true),:xaxis => KW(:autorange => true)));
	annotate!((real(constellationQ132[:]),imag(constellationQ132[:]) .+ .06,alphabetQ132[:],10),title="4-QAM");
end

# ╔═╡ d087448c-8e8d-11eb-20f2-a77402b7fde6
html"""
<h4>Question 2</h4>
<p>We have the exact error probability for a rectangular M-QAM is given by <a href="#ref-proakis">[Proakis]</a> at Equation 4.3-30, just noting that $\mathcal{E}_\text{avg} = \log_{10}M\mathcal{E}_\text{bavg}$

$$\begin{aligned}
P_{e, M-\mathrm{QAM}}=& 4\left(1-\frac{1}{\sqrt{M}}\right) Q\left(\sqrt{\frac{3 }{M-1} \frac{\mathcal{E}_{\mathrm{avg}}}{N_{0}}}\right) \\
& \times\left(1-\left(1-\frac{1}{\sqrt{M}}\right) Q\left(\sqrt{\frac{3 }{M-1} \frac{\mathcal{E}_{\mathrm{avg}}}{N_{0}}}\right)\right)
\end{aligned}.$$

If we consider the $\text{SNR}=\mathcal{E}_{\mathrm{avg}}/N_{0}$ per bit we obtain the curves below. The <a href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/BERSER.jl" target="_blank">implementation</a> at the function <code>TER</code> gives the <strong>Theoretical Error Rates</strong> in function of the <strong>Complementary Error Function</strong>, using the relation $$ Q(x) = \frac{1}{2}\operatorname{erfc} \left(\frac{x}{\sqrt{2}} \right) $$ 

Then we plot in the interval <code>SNRdB = 0:1:20</code>.</p>
"""

# ╔═╡ a4585958-8f41-11eb-374d-056505794fdb
begin
	plotly();
	SNRdB = 0:1:20;
	SER4QAM,BER4QAM = TER("QAM",4,SNRdB);
	SER16QAM,BER16QAM = TER("QAM",16,SNRdB);
	SER64QAM,BER64QAM = TER("QAM",64,SNRdB);
	p = plot(SNRdB,[SER64QAM,SER16QAM,SER4QAM],label=["64-QAM" "16-QAM" "4-QAM"],legend=(.3,.3),lw=2,yscale=:log10,ylim=(1e-6,1),xlim=(0,20),yminorticks=5);
	xticks!(collect(0:2:20),title="Theoretical SER",xlabel="SNR(dB)",ylabel="SER");
end

# ╔═╡ ebaa895c-8f7a-11eb-3a19-d3b8330dbc01
md"""
#### Question 3

The **Symbol Error Rate** (SER), as the name suggests, it's the rate between the number of wrong symbols and the number of sent symbols. Similarly is the **Bit Error Rate**, but looking each bit individually. To simulate the wrong bits, we need to explain how the process of sending a symbol happens in this work.

###### 1. Generating the message
First of all we need to generate the message to be sent using one of our modulations. We generate a sequence of bits `bitseq` in the line below.
"""

# ╔═╡ c4fa94ce-9007-11eb-166d-bd7f97682161
bitseq = BitArray(rand((0,1),1,264000));

# ╔═╡ 38f8742e-9018-11eb-0b6e-dfced0559513
html"""
<h6>2. Generating the symbol constelation</h6>
<p>
After we need to indicate where are the symbols of our modulation are. We can do this using one of the modulation functions <code>MQAM</code> or <code>MPSK</code> (<i>the code for PSK will be presented below when we treat about it</i>).

This information can be passed directly to <code>bitModulation</code> function implemented <a href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/bitModulation.jl" target="_blank">here</a>. The function creates a <code>dictionary</code> that assign each symbol to its respective code word in bits. Then we slice the initial <code>bitseq</code> in <code>chunks</code> of the size <code>M</code> of the modulation. Finally we assign the respective symbol to each chunck.</p>
"""

# ╔═╡ 6674fe84-901a-11eb-11a5-7fcb2e63f29b
md"""
	signal4QAM = bitModulation(bitseq,MQAM(4));

As shown above, `bitModulation` modulates `bitseq` with a 4-QAM, i.e `MQAM(4)`. In the code below, we do this for 4, 16 and 64 orders of the QAM"
"""

# ╔═╡ a447330c-9007-11eb-2467-27570bf1c6f2
begin
	signal4QAM = bitModulation(bitseq,MQAM(4));
	signal16QAM = bitModulation(bitseq,MQAM(16));
	signal64QAM = bitModulation(bitseq,MQAM(64));
end

# ╔═╡ b6ad2b5c-9019-11eb-0f45-0707b994a43e
md"""
We can see an example of some symbols for the `signal16QAM` variable assigned above.
"""

# ╔═╡ edefcb92-9019-11eb-0340-a513618513b3
scatter(signal16QAM[1:100],title="16-QAM Symbols",legend=:false)

# ╔═╡ 84219e1e-9020-11eb-31fe-d526c59b81c6
html"""
<h6>3. Simulating a channel</h6>
<p>The next step is sent the <code>signal16QAM</code> through a <code>AWGN</code> channel, implemented <a href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/AWGN.jl" target="_blanck">here</a>. The noise added by the channel has <b>standard deviation</b> $\sigma = \sqrt{N_0/2}$, which can be passed in function of <code>√avgEnergy</code> which is the <b>average symbol energy</b> of the sent signal <code>σ = √(avgEnergy/SNR)</code>. The <code>SNR</code>(Signal Noise Ratio). We <b>rescale</b> the signal once sent to recover its original constelation position <a href="#ref-dsplog">[DSPLOG]</a>, multiplying it by <code>√avgEnergy</code>. The sinal after the channel is plotted below with a <code>SNR</code> of 20 dB.</p>
"""

# ╔═╡ 3012015e-9022-11eb-162b-714fbae56ef1
scatter(AWGN(signal16QAM[1:100],10.0.^(20/10)),title="16-QAM Symbols",legend=:false)

# ╔═╡ b59bed5c-9023-11eb-3779-63f9c1cb8f9c
html"""
<h6>4. Demodulating the received signal</h6>
<p>Once received, we must pass the signal through a <b>classifier</b> which decides to which symbol in the original constellation the received symbol belongs. As the symbols in the original constellation has the <b>same probability</b>, the <b> optimal threshold</b> is at the <b>half of the distance between neighbor symbols</b>. This is <a target="_blank" href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/signalClassifier.jl">implemented</a> in the <code>absClassifier</code>, where we evaluate the <b>absolute distance</b> in</p>

<pre><code>abs.(transpose(constellation[:]) .- signal[:])</code></pre>

<p>between every sent symbol and each symbol in the original constellation. The <code>findmin</code> function finds, at each sent symbol, the <b>closer constellation symbol</b>. Wrapping up, this distance is <b>always smaller</b> at the side of the threshold for which symbol we have to decide for. This is similar to "slice" the space where the symbols are, according to where the threshold lays on, which is exactly where the distance is equal for the both neighbor original constellation symbols.</p>

<p>Finally, the decided symbols are assigned to its respective code word, at the <code>bitDemodulation</code> dunction implemented <a target="_blank" href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/bitDemodulation.jl">here</a>.</p>
"""

# ╔═╡ a54b1416-9026-11eb-3560-a3c9bd59129d
html"""
<h6>5. Counting the errors</h6>
<p>In this last part, we count the bits which the classifier decided wrong. Using the function <code>BER</code>, implemented <a href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/BERSER.jl" target="_blank">here</a>, we compare the sent message <code>signal4QAM</code> with the one after the channel in <code>sum(bitseq .!= rbitseq)</code>. For reasons of implementation, all the previous steps <b>3</b> and <b>4</b> were applied just at this point, as can be saw in the code.</p>

<p> This function can implement too a <b>Monte Carlo</b> simulation. We pass the signal <code>iter</code> times through the <code>AWGN</code> channel. This may be necessary in simulations fow AWGN channels with SNR greater than 15 dB, as experienced at this work. This due the fact that, for detect a wrong signal, it must lay on the wrong side of the threshold. But if the noise is not strong enough, the probability of it pushes the symbol beyond the threshold is too small, then we need to sample several time in order to obtain a number closer to the theoretical result.</p>

<p>Below, the <code>Simu*</code> variables receive the results of the <code>BER</code> and <code>SER</code> functions.
"""

# ╔═╡ 2a40e80c-900a-11eb-29a8-31ca61a382f4
begin
		SimuBER4QAM = map(snr-> BER(signal4QAM,snr,MQAM(4),iter=10),SNRdB);
		SimuBER16QAM = map(snr-> BER(signal16QAM,snr,MQAM(16)),SNRdB);
		SimuBER64QAM = map(snr-> BER(signal64QAM,snr,MQAM(64)),SNRdB);
		
		SimuSER4QAM = map(snr-> SER(signal4QAM,snr,MQAM(4),iter=10),SNRdB);
		SimuSER16QAM = map(snr-> SER(signal16QAM,snr,MQAM(16)),SNRdB);
		SimuSER64QAM = map(snr-> SER(signal64QAM,snr,MQAM(64)),SNRdB);
end

# ╔═╡ 0afb3900-9007-11eb-2517-d1d6eb9788f2
begin
	pQ3BER = plot(SNRdB,[SimuBER4QAM SimuBER16QAM SimuBER64QAM]);
	pQ3SER = plot(SNRdB,[SimuSER4QAM,SimuSER16QAM,SimuSER64QAM]);
	plot(pQ3SER,yaxis=:log10,ylim=(1e-5,1),title="Simulated SER",xlabel="SNR(dB)",ylabel=["SER" "BER"],label=["4-QAM" "16-QAM" "64-QAM"],mark=:o,extra_plot_kwargs = KW(:yaxis => KW(:autorange => true),:xaxis => KW(:autorange => true)));
	xticks!(collect(0:2:20));
end

# ╔═╡ cbb6ab50-900e-11eb-09a8-d5b2292bb292
begin
	plot(pQ3BER,yaxis=:log10,ylim=(1e-5,1),xlabel="SNR(dB)",ylabel=["BER" "BER"],title="Simulated BER",label=["4-QAM" "16-QAM" "64-QAM"],mark=:o,extra_plot_kwargs = KW(:yaxis => KW(:autorange => true),:xaxis => KW(:autorange => true)))
	xticks!(collect(0:2:20));
	
	
end

# ╔═╡ 1cb43388-902a-11eb-0972-8d63b7535172
md"""
#### Question 4.1

###### 1. Constellation Average Energy

The Average Constellation Energy for the PSK can be obtained in a similar way of the QAM. We consider the same basis functions, just **adding a phase shifting term** [[Proakis]](#ref-proakis), obtaining the bidimensional form, 

$\mathbf{s}_m=\left(\sqrt{\frac{\mathcal{E}_{g}}{2}} \cos \left(\frac{2 \pi}{M}(m-1)\right), \sqrt{\frac{\mathcal{E}_{g}}{2}} \sin \left(\frac{2 \pi}{M}(m-1)\right)\right), \quad m=1,2,\dots,M.$

Note the norm of the symbol does **not** change, always laying in a circle of square radius equals the half of the pulse energy, then 

$$\mathcal{E}_\text{avg} = \frac{1}{M} \cdot \sum_{m=1}^M \frac{\mathcal{E}_g}{2} = \frac{1}{2}\mathcal{E}_g.$$

"""

# ╔═╡ 7300fd16-9032-11eb-19b2-096883eb3511
md"""
###### 2. Minimum distance

By the bidimensional form given above, we have that the minumum distance between signal points is

$\begin{aligned}
d_{m n} &=\sqrt{\left\|\mathbf{s}_{m}-\mathbf{s}_{n}\right\|^{2}} \\
&=\sqrt{\mathcal{E}_{g}\left[1-\cos \left(\frac{2 \pi}{M}(m-n)\right)\right]}
\end{aligned}$

If $n$ is the neighbor, then $|m-n|=1$. Knowing that $\cos (2 \theta)=1-2 \sin ^{2} \theta$, we obtaing

$d_{\min }=\sqrt{\mathcal{E}_{g}\left(1-\cos \frac{2 \pi}{M}\right)}=\sqrt{2 \mathcal{E}_{g} \sin ^{2} \frac{\pi}{M}}.$

"""

# ╔═╡ 734fd4da-9033-11eb-1d25-ebe0ecfdf570
html"""
<h6> 3. Modulator and demodulator</h6>

<p>Both modulator and demodulator are implemented in a similar way of the QAM, just using the symbols generated by the <code>PSK</code> function implemented <a target="_blank" href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/PSK.jl">here</a>. It uses the fact of the bidimensional form above as in

<pre><code>phase = 2pi/M .* (m.-1);
constellation = √(.5energy) * ( cos.(phase) + 1im .* sin.(phase));</code></pre>
"""

# ╔═╡ f8c2c728-902c-11eb-3d65-9300fe828483
begin
 	m,alphabetQ4,constellationQ4 = MPSK(4,1,unitAveragePower=false);
	scatter(real(constellationQ4[:]),imag(constellationQ4[:]),legend=false,extra_plot_kwargs = KW(:yaxis => KW(:autorange => true),:xaxis => KW(:autorange => true)));
	annotate!((real(constellationQ4[:]),imag(constellationQ4[:]) .+ .06,alphabetQ4[:],10));
end

# ╔═╡ 08508e66-9035-11eb-3b5e-41f34f2c4a50
begin
 	m42,alphabetQ42,constellationQ42 = MPSK(8,1,unitAveragePower=false);
	scatter(real(constellationQ42[:]),imag(constellationQ42[:]),legend=false,extra_plot_kwargs = KW(:yaxis => KW(:autorange => true),:xaxis => KW(:autorange => true)));
	annotate!((real(constellationQ42[:]),imag(constellationQ42[:]) .+ .06,alphabetQ42[:],10));
end

# ╔═╡ c838a1a0-9035-11eb-3efb-ff796ddf6d25
html"""
<h4>Question 4.2-3</h4>
<p>As in<b>Question 2</b> the <a href="https://raw.githubusercontent.com/filippfarias/DigCommSys/master/src/BERSER.jl" target="_blank">implementation</a> at the function <code>TER</code> gives the <strong>Theoretical Error Rates</strong> using the expression <a href="#ref-proakis">[Proakis]</a> at Equation 4.3-34, just noting that $\mathcal{E}_\text{avg} = \log_{10}M\mathcal{E}_\text{bavg}$

$$P_{M} \approx 2 Q\left(\sqrt{\left(2 \log _{2} M\right) \sin ^{2}\left(\frac{\pi}{M}\right) \frac{\mathcal{E}_{b}}{N_{0}}}\right).$$

Then we plot in the interval <code>SNRdB = 0:1:20</code>.</p>
"""

# ╔═╡ 8624d396-9036-11eb-3bd0-216cafe61255
begin
	SER4PSK,BER4PSK = TER("PSK",4,SNRdB);
	SER8PSK,BER8PSK = TER("PSK",8,SNRdB);
	p42 = plot(SNRdB,[SER4PSK,SER8PSK],label=["4-PSK" "8-PSK"],legend=(.3,.3),lw=2,yscale=:log10,ylim=(1e-6,1),xlim=(0,20));
	xticks!(collect(0:2:20),title="Theoretical SER",xlabel="SNR(dB)",ylabel="SER");
end

# ╔═╡ 7550330e-9010-11eb-1cba-695d74093019
begin
	signal4PSK = bitModulation(bitseq,MPSK(4));
	signal8PSK = bitModulation(bitseq,MPSK(8));
end

# ╔═╡ 726600fe-900f-11eb-2c5a-23c6e9d2da14
begin
		SimuBER4PSK = map(snr-> BER(signal4PSK,snr,MPSK(4)),SNRdB);
		SimuBER8PSK = map(snr-> BER(signal8PSK,snr,MPSK(8)),SNRdB);
		
		SimuSER4PSK = map(snr-> SER(signal4PSK,snr,MPSK(4)),SNRdB);
		SimuSER8PSK = map(snr-> SER(signal8PSK,snr,MPSK(8)),SNRdB);
end

# ╔═╡ da42142e-900f-11eb-3b2e-450009059ed3
begin
	pQ41BER = plot(SNRdB,[SimuBER4PSK SimuBER8PSK]);
	pQ41SER = plot(SNRdB,[SimuSER4PSK,SimuSER8PSK]);
	plot(pQ41SER,yaxis=:log10,ylim=(1e-5,1),title="Simulated SER",xlabel="SNR(dB)",ylabel=["SER" "BER"],label=["4-PSK" "8-PSK"],mark=:o,extra_plot_kwargs = KW(:yaxis => KW(:autorange => true),:xaxis => KW(:autorange => true)));
end

# ╔═╡ dac483be-900f-11eb-02e1-49b0840b3798
plot(pQ41BER,yaxis=:log10,ylim=(1e-5,1),xlabel="SNR(dB)",ylabel=["BER" "BER"],title="Simulated BER",label=["4-PSK" "8-PSK"],mark=:o,extra_plot_kwargs = KW(:yaxis => KW(:autorange => true),:xaxis => KW(:autorange => true)))

# ╔═╡ 9923382e-9037-11eb-26b7-75e85a4171c5
md"""
#### Question 5.1
With ball marks, the simulated.
"""

# ╔═╡ ae9ec3ee-9037-11eb-0fc4-3329d863e33d
begin
	plot(SNRdB,[SER64QAM,SER16QAM,SER4QAM,SER4PSK,SER8PSK],label=["64-QAM" "16-QAM" "4-QAM" "4-PSK" "8-PSK"],legend=(.3,.3),lw=2,yscale=:log10,ylim=(1e-6,1),xlim=(0,20),yminorticks=5);
	plot!(SNRdB,[SimuSER64QAM,SimuSER16QAM,SimuSER4QAM,SimuSER4PSK,SimuSER8PSK],label=["64-QAM" "16-QAM" "4-QAM" "4-PSK" "8-PSK"],legend=(.3,.3),lw=2,yscale=:log10,ylim=(1e-6,1),xlim=(0,20),mark=:o,markercolor = :transparent);
	xticks!(collect(0:2:20),title="SER",xlabel="SNR(dB)",ylabel="SER");
end

# ╔═╡ 0db87482-9039-11eb-08bc-2593ecb43c37
md"""
#### Question 5.2

With ball marks, the simulated.
"""

# ╔═╡ d79805ca-9038-11eb-2479-8de85915601b
begin
	plot(SNRdB,[BER64QAM,BER16QAM,BER4QAM,BER4PSK,BER8PSK],label=["64-QAM" "16-QAM" "4-QAM" "4-PSK" "8-PSK"],legend=(.3,.3),lw=2,yscale=:log10,ylim=(1e-6,1),xlim=(0,20),yminorticks=5);
	plot!(SNRdB,[SimuBER64QAM,SimuBER16QAM,SimuBER4QAM,SimuBER4PSK,SimuBER8PSK],label=["64-QAM" "16-QAM" "4-QAM" "4-PSK" "8-PSK"],legend=(.3,.3),lw=2,yscale=:log10,ylim=(1e-6,1),xlim=(0,20),mark=:o,markercolor = :transparent);
	xticks!(collect(0:2:20),title="BER",xlabel="SNR(dB)",ylabel="BER");
end

# ╔═╡ 10acaf46-9039-11eb-0943-e93ab2b8ca9a
md"""
#### Question 5.3
In the later Plots we see that 4-QAM and 4-PSK have the **similar efficiencies**, what could be expecter given that one is just the phase shift of another. It is important to note the **decrease of the error rate** as the order of the modulation **rises**, but it requires a larger **update rate** when transmitting given that we send more bits per symbol, implying in a **larger frequency band**. We see to the increase of the SNR, i.e. the energy necessary to overcome the effect of the noise given that, when normalized to 1, the symbols become closer.
"""

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
<p><a name="ref-dsplog">[DSPLOG]</a> <a target="_blank" href="http://www.dsplog.com/2012/01/01/symbol-error-rate-16qam-64qam-256qam/">Symbol Error rate for QAM (16, 64, 256,.., M-QAM)</a>. DSPLOG - Signal Processing for Communication.</p>
<p><a name="ref-ieee">[IEEE]</a> <a target="_blank" href="https://en.wikipedia.org/wiki/Gray_code#Constructing_an_n-bit_Gray_code">Bit error probability of M-ary quadrature amplitude modulation</a>. Dongweon Yoon, Kyongkuk Cho, Jinsock Lee. </p>
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
Note the plot below is **not normalised** by the average energy, but the normalization is required when the signal is sent.""");
end

# ╔═╡ 2a7fc380-9032-11eb-3eab-9130582c4c48
begin 
	M2_ = [4,8];
	e2_avg = 1 ./ 2;
	correct("Answer for Question 4.1",md"""
For $M=\{4,8\}$, considering $\mathcal{E}_g=1$, the average energy is $ $[e2_avg]. $ for both.""");
end

# ╔═╡ 7613d314-9034-11eb-332d-17ffdf9712bb
begin 
	correct("Answer for Question 4.2",md"""
If $\sqrt{\mathcal{E}_g}=1$, then for $M=\{4,8\}$ is equal to $\sqrt{2}$ and $1$, respectively.""");
end

# ╔═╡ Cell order:
# ╟─8df467aa-8bd9-11eb-088e-218f062da936
# ╠═25bc3a22-8be4-11eb-1414-0b8e408a7b54
# ╠═2780cb8e-901c-11eb-3bf6-6149f19d97c3
# ╟─f1835c24-8bdc-11eb-2e3e-ffaca518a99c
# ╟─57549df8-8c3f-11eb-16a5-c33fe4560e41
# ╟─9dbdeb40-8c24-11eb-1569-6f010c7f5ed5
# ╟─98e4d1ba-8c45-11eb-2e90-21ecfb046773
# ╠═7dfecccc-8c5c-11eb-01e5-1f9b7f06dbc8
# ╟─df072ef6-8c63-11eb-0064-8328864e88cf
# ╟─8036d610-8c61-11eb-2cef-43773fdbf97c
# ╟─53d7bd68-8f3a-11eb-251c-d9b4a34b8a44
# ╟─bd6f75c6-8c63-11eb-1541-f17f0f1803bb
# ╟─45b663b0-8c67-11eb-1302-df7ec6e9452d
# ╟─684ef8fc-8c69-11eb-0947-771929a7f355
# ╟─4d3d96b6-8ce0-11eb-21c2-ebfbda60351f
# ╠═8d99049a-8c6b-11eb-3118-51c16e0927bf
# ╟─78f40b68-8c6c-11eb-2925-15974d7dccb5
# ╟─671b4d78-9035-11eb-0a9e-a3abea476a77
# ╟─68a8d124-9035-11eb-0d05-97e7e4454d10
# ╟─d087448c-8e8d-11eb-20f2-a77402b7fde6
# ╠═a4585958-8f41-11eb-374d-056505794fdb
# ╟─ebaa895c-8f7a-11eb-3a19-d3b8330dbc01
# ╠═c4fa94ce-9007-11eb-166d-bd7f97682161
# ╟─38f8742e-9018-11eb-0b6e-dfced0559513
# ╟─6674fe84-901a-11eb-11a5-7fcb2e63f29b
# ╠═a447330c-9007-11eb-2467-27570bf1c6f2
# ╟─b6ad2b5c-9019-11eb-0f45-0707b994a43e
# ╠═edefcb92-9019-11eb-0340-a513618513b3
# ╟─84219e1e-9020-11eb-31fe-d526c59b81c6
# ╠═3012015e-9022-11eb-162b-714fbae56ef1
# ╟─b59bed5c-9023-11eb-3779-63f9c1cb8f9c
# ╟─a54b1416-9026-11eb-3560-a3c9bd59129d
# ╠═2a40e80c-900a-11eb-29a8-31ca61a382f4
# ╟─0afb3900-9007-11eb-2517-d1d6eb9788f2
# ╟─cbb6ab50-900e-11eb-09a8-d5b2292bb292
# ╟─1cb43388-902a-11eb-0972-8d63b7535172
# ╟─2a7fc380-9032-11eb-3eab-9130582c4c48
# ╟─7300fd16-9032-11eb-19b2-096883eb3511
# ╟─7613d314-9034-11eb-332d-17ffdf9712bb
# ╟─734fd4da-9033-11eb-1d25-ebe0ecfdf570
# ╟─f8c2c728-902c-11eb-3d65-9300fe828483
# ╟─08508e66-9035-11eb-3b5e-41f34f2c4a50
# ╟─c838a1a0-9035-11eb-3efb-ff796ddf6d25
# ╠═8624d396-9036-11eb-3bd0-216cafe61255
# ╠═7550330e-9010-11eb-1cba-695d74093019
# ╠═726600fe-900f-11eb-2c5a-23c6e9d2da14
# ╠═da42142e-900f-11eb-3b2e-450009059ed3
# ╠═dac483be-900f-11eb-02e1-49b0840b3798
# ╟─9923382e-9037-11eb-26b7-75e85a4171c5
# ╠═ae9ec3ee-9037-11eb-0fc4-3329d863e33d
# ╟─0db87482-9039-11eb-08bc-2593ecb43c37
# ╠═d79805ca-9038-11eb-2479-8de85915601b
# ╟─10acaf46-9039-11eb-0943-e93ab2b8ca9a
# ╟─0d70be04-8bf6-11eb-0638-a35109c2d36c
# ╟─21953290-8bf7-11eb-025f-87f8477b88c2
# ╟─8e37ac06-8c5b-11eb-3ca4-07b42d5b841a
# ╟─bf54e3d6-8c5b-11eb-1cd9-2b7bc38fc404
# ╟─bf554652-8c5b-11eb-30f1-b9352a6f77cf
# ╟─bf5fc848-8c5b-11eb-16e1-b7c9b467c87a
