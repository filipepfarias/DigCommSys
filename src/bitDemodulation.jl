"""
---
bitDemodulation(signal,M)
---
"""
function bitDemodulation(signal,M)
    m = Int(log2(M))
    N = length(signal)*m
    if mod(N,M) == 0
        (alphabet, constellation) = MQAM(M)

        d = 2; # To define

        dictionary = Dict( 
            constellation[i] => alphabet[i] for i=1:length(constellation)
            )
        
        (_,ids) = findmin(
            abs.(transpose(constellation[:]) .- signal[:]),
            dims=2)

        chunks = map(p -> dictionary[constellation[p[2]]],ids)
        
        bitsequence = BitArray(transpose(
            BitArray(map(
                bit -> parse(Int,bit), 
                collect(prod(chunks))
                ))
        ))

        return bitsequence
    else
        error("The bit sequence is not multiple of M.")
    end
end
