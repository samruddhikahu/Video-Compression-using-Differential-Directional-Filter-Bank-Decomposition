function MVbits = MVD_enc(MV)
p = size(MV,2);
q = size(MV,3);
y = 1; x = 1;
MVbits{1} = '';
MVbits{2} = '';
i = 1;
while(i<=(p*q))
    clear bits;
    bits = enc_golomb(MV(1,y,x),1);
    MVbits{1} = [MVbits{1} bits];
    clear bits;
    bits = enc_golomb(MV(2,y,x),1);
    MVbits{2} = [MVbits{2} bits];
    if(mod(i,p)==0)
        x = x + 1;
        y = 1;
    else
        y = y + 1;
    end
    i = i + 1;
end

end