function MV_rec = MVD_dec(MVbits,mv_fr,nfr)
MV_rec = zeros(2,mv_fr,nfr);
y = 1; x = 1;
idx1 = 1; idx2 = 1;
i = 1;
while(i<=(mv_fr*nfr))
    [MV_rec(1,y,x),idx1] = dec_golomb(idx1,MVbits{1},1);
    [MV_rec(2,y,x),idx2] = dec_golomb(idx2,MVbits{2},1);
    if(mod(i,mv_fr)==0)
        x = x + 1;
        y = 1;
    else
        y = y + 1;
    end
    i = i + 1;
end
end