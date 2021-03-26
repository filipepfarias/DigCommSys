"""
---
BER
---
"""
function BER(M,signal,SNR,iter=1)
    i = 1;
    Eber = ones(iter,1);

    while i<=iter
        tsignal = AWGN(signal,SNR);
        rbitseq = bitDemodulation(tsignal,M);
        bitseq = bitDemodulation(signal,M);

        Eber[i] = sum(bitseq .== rbitseq)/length(bitseq);

        i += 1;
    end
    return Eber
end

"""
---
SER
---
"""
function SER(M,signal,SNR,iter=1)
    i = 1;
    Eser = ones(iter,1);

    while i<=iter
        tsignal = AWGN(signal,SNR);
        rbitseq = bitDemodulation(tsignal,M);
        dsignal = bitModulation(rbitseq,M);

        Eser[i] = sum(signal .== dsignal)/length(signal);

        i += 1;
    end
    return Eser
end

"""
---
MCER = Monte Carlo Error Rate
---
"""
function MCER(ER)
    return cumsum(ER,dims=1)./(1:length(ER));
end