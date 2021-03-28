"""
---
PSK Modulation
---
"""
function MPSK(M,energy=2;unitAveragePower=true)
    if M == 8 || M == 4
        m = 1:M;
        phase = 2pi/M .* (m.-1);
        avgEnergy = .5energy;

        constellation = √(.5energy) * ( cos.(phase) + 1im .* sin.(phase));

        unitAveragePower ? constellation = constellation ./ √avgEnergy : nothing;

        alphabet = GrayCode(Int(log2(M)));
        alphabet = [bitstring(alphabet[i,:]) for i=1:size(alphabet)[1]];

        return M, alphabet, constellation
    else
        error("Select M equals 4 or 8.")
    end
end