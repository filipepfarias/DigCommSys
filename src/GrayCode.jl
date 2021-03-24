"""
GrayCode(N)
---
Generates Gray Code lookup table for N bits.
"""
function GrayCode(N)
	seq = BitArray([0, 1])
	i = 1
	while i<N
		Z = BitArray(zeros(size(seq)[1]))
		seq = vcat(hcat(Z,seq),hcat(.~Z,seq[end:-1:1,:]))
		i = i+1
	end
	return seq
end