"""
---
AWGN Channel
---
"""
function AWGN(signal,SNR)
    N = length(signal);

    avgEnergy = sum( abs.(signal[:]).^2 )/N;

    σ = √(avgEnergy/SNR);

    tsignal = signal.*√avgEnergy + σ/√2 * (randn(N,1) + 1im*randn(N,1));
    
    return tsignal
end