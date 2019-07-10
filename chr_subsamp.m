function Cs = chr_subsamp(C)
[m, n] = size(C);
Cs = zeros(m/2,n/2);
k = 1;
for i = 1:2:m
    l = 1;
    for j = 1:2:n
        G = C(i:i+1,j:j+1);
        Avg = (sum(sum(G)))/4;
        Cs(k,l) = Avg;
        l = l + 1;
    end
    k = k + 1;
end
end