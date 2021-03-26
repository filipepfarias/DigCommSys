"""
---
bitModulation(bitsequence,M)
---
"""
function bitModulation(bitsequence::BitArray,(M,alphabet, constellation))
    N = length(bitsequence)
    m = Int(log2(M))
    if mod(N,M) == 0
        dictionary = Dict( alphabet[i] => constellation[i] for i=1:length(constellation))

        chunks = transpose(reshape(bitsequence,(m,floor(Int,N/m))))

        signal = [dictionary[bitstring(chunks[i,:])] for i=1:size(chunks)[1]]
        return signal
    else
        error("The bit sequence is not multiple of M.")
    end
end