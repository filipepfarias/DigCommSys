"""
---
QAM Modulation
---
"""
function MQAM(M,energy=2,unitAveragePower=false)
    if M==sqrt(M)^2 && 0==mod(M,2)
        m = Int(√M);
        d = sqrt(2*energy); #Digital Communications-Proakis (Table 3.2-1)
        avgEnergy = (M-1)/3*energy;

        I = ((2*(1:m) .- 1 .- m)*d/2) .+0*1im;
        Q = I.*1im
        constellation = (Q.+transpose(I));
        
        unitAveragePower ? constellation = constellation ./ √avgEnergy : nothing;

        alphabet = GrayCode(Int(log2(M)))
        alphabet = [bitstring(alphabet[i,:]) for i=1:size(alphabet)[1]]
        alphabet = reshape(alphabet,(m,m))

        for i=2:2:m
            alphabet[:,i] = alphabet[end:-1:1,i]
        end

        return M, alphabet, constellation
    else
        error("M is neither square nor even.")
    end
end