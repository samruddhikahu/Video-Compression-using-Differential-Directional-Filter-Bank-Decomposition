function Sband = rearrange1(coeffs, lev)

if(lev==1)
    Sband = [coeffs{1}; coeffs{2}];
else
    Sb1 = []; Sb2 = []; Sb3 = []; Sb4 = [];
    xx = 2^(lev-2);
    for k = 1:xx
        Sb1 = [Sb1; coeffs{k}];
        Sb2 = [Sb2; coeffs{xx + k}];
        Sb3 = [Sb3 coeffs{2*xx + k}];
        Sb4 = [Sb4 coeffs{3*xx + k}];
    end
    Sband = [Sb1 Sb2;
        Sb3 Sb4];
end

end