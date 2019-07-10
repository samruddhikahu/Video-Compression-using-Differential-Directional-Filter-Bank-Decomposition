function coeffs = rearrange_back(Sbands, lev)
[m,n] = size(Sbands);

if(lev==1)
    coeffs{1} = Sbands(1:(m/2),1:n);
    coeffs{2} = Sbands((m/2)+1:m,1:n);
else
    xx = 2^(lev-2);
    Sb1 = Sbands(1:(m/2),1:(n/2));
    Sb2 = Sbands(1:(m/2),(n/2)+1:n);
    Sb3 = Sbands((m/2)+1:m,1:(n/2));
    Sb4 = Sbands((m/2)+1:m,(n/2)+1:n);
    a = 2*xx;
    for k = 1:xx
        coeffs{k} = Sb1(1:(m/a),:);
        Sb1(1:(m/a),:) = [];
        coeffs{xx + k} = Sb2(1:(m/a),:);
        Sb2(1:(m/a),:) = [];
        coeffs{2*xx + k} = Sb3(:,1:(n/a));
        Sb3(:,1:(n/a)) = [];
        coeffs{3*xx + k} = Sb4(:,1:(n/a));
        Sb4(:,1:(n/a)) = [];
    end
end

end