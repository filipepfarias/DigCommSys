"""
---
PSK Modulation
---
"""
function MPSK(M,energy=2)
    if M == 8 || M == 4
        m = 1:M;
        phase = 2pi/M .* (m.-1);

        constellation = √(energy/2) * ( cos.(phase) + 1im .* sin.(phase));

        alphabet = GrayCode(Int(log2(M)));
        alphabet = [bitstring(alphabet[i,:]) for i=1:size(alphabet)[1]];

        return M, alphabet, constellation
    else
        error("Select M equals 4 or 8.")
    end
end