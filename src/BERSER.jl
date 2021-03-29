"""
---
BER
---
"""
function BER(signal,SNRdB,(M,alphabet,constellation);iter=1,saveResults=false)
    i = 1;
    Eber = ones(iter,1);
    SNR = 10 .^(SNRdB/10);

    while i<=iter
        tsignal = AWGN(signal,SNR);
        rbitseq = bitDemodulation(tsignal,(M,alphabet,constellation));
        bitseq = bitDemodulation(signal,(M,alphabet,constellation));

        if iter > 1
            Eber[i] = BigFloat(sum(bitseq .!= rbitseq)/length(bitseq));
        else
            Eber = sum(bitseq .!= rbitseq)/length(bitseq);
        end
        i += 1;
    end
    saveResults ? Eber : Eber = sum(Eber)/length(Eber);
    return Eber
end

"""
---
SER
---
"""
function SER(signal,SNRdB,(M,alphabet,constellation);iter=1,saveResults=false)
    i = 1;
    Eser = ones(iter,1);
    SNR = 10 .^(SNRdB/10);

    while i<=iter
        tsignal = AWGN(signal,SNR);
        rbitseq = bitDemodulation(tsignal,(M,alphabet,constellation));
        dsignal = bitModulation(rbitseq,(M,alphabet,constellation));

        if iter > 1
            Eser[i] = sum(signal .!= dsignal)/length(signal);
        else
            Eser = sum(signal .!= dsignal)/length(signal);
        end
        i += 1;
    end
    saveResults ? Eser : Eser = sum(Eser)/length(Eser);
    return Eser
end

"""
---
MCER - Monte Carlo Error Rate
---
"""
function MCER(ER)
    return cumsum(ER,dims=1)./(1:length(ER));
end

"""
---
Theorical Error Rates
---
"""
function TER(mod,M,SNRdB)
    SNR = 10 .^(SNRdB/10);
    b = log2(M);
    
    if mod == "QAM"
        q = 4*(√M-1)/√M*Q(sqrt.(3/(M-1)*SNR));
        SER = (1 .- q/4) .* q;
        if (M == 4)
            BER = Q(sqrt.(2*SNR/b));
        elseif (M == 16)
            BER = 3/4*Q(sqrt.(4/5*SNR/b)) + 1/2*Q(3*sqrt.(4/5*SNR/b)) - 1/4*Q(5*sqrt.(4/5*SNR/b));
        elseif (M == 64)
            BER = 7/12*Q(sqrt.(2/7*SNR/b)) + 1/2*Q(3*sqrt.(2/7*SNR/b)) - 1/12*Q(5*sqrt.(2/7*SNR/b)) + 1/12*Q(9*sqrt.(2/7*SNR/b)) - 1/12*Q(13*sqrt.(2/7*SNR/b));
        end
    elseif mod == "PSK"
        SER = 2Q(sqrt.( 2 * SNR ) * sin(pi/M));
        BER = SER/b;
    else
        error("Not valid modulation.")
    end
    return SER,BER
end

"""
---
Q Function
---
"""
function Q(x)
    return .5* SpecialFunctions.erfc.(x/√2)
end