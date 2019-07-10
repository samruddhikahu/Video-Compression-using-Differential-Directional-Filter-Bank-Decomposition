function C = chr_resamp(Cs)
[m1, n1] = size(Cs);
k = 1;
for i = 1:2:(m1*2)
    l = 1;
    for j = 1:2:(n1*2)
        C(i:i+1,j:j+1) = Cs(k,l);
        l = l + 1;
    end
    k = k + 1;
end
end