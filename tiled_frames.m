function Frs = tiled_frames(Fr1,Fr2,Fr3,TileSize)
[m1,n1,o,p] = size(Fr1);
[m2,n2] = size(Fr2(:,:,1,1));
[m3,n3] = size(Fr3(:,:,1,1));

mt = round(TileSize(1,1));
nt = round(TileSize(1,2));

if(mt~=TileSize(1,1))
    x1 = m1/TileSize(1,1);
    Fr1((m1+1):(mt*x1),:,1,:) = 0;
    
    x2 = m2/TileSize(1,1);
    Fr2((m2+1):(mt*x2),:,1,:) = 0;
    Fr3((m3+1):(mt*x2),:,1,:) = 0;
end
if(nt~=TileSize(1,2))
    y1 = n1/TileSize(1,2);
    Fr1(:,(n1+1):(nt*y1),1,:) = 0;
    
    y2 = n2/TileSize(1,2);
    Fr2(:,(n2+1):(nt*y2),1,:) = 0;
    Fr3(:,(n3+1):(nt*y2),1,:) = 0;
end
[m1,n1,o,p] = size(Fr1);
[m2,n2] = size(Fr2(:,:,1,1));
[m3,n3] = size(Fr3(:,:,1,1));

nom1 = m1/mt; non1 = n1/nt;
nom2 = m2/mt; non2 = n2/nt;
nom3 = m3/mt; non3 = n3/nt;

k = 1;
for ii = 1:p
    for i = 1:nom1
        for j = 1:non1
            a1 = (i-1)*mt + 1; b1 = i*mt;
            a2 = (j-1)*nt + 1; b2 = j*nt;
            Frs(:,:,1,k) = Fr1(a1:b1,a2:b2,1,ii);
            k = k + 1;
        end
    end
    if((nom2>=1)||(non2>=1))
        for i = 1:nom2
            for j = 1:non2
                a1 = (i-1)*mt + 1; b1 = i*mt;
                a2 = (j-1)*nt + 1; b2 = j*nt;
                Frs(:,:,1,k) = Fr2(a1:b1,a2:b2,1,ii);
                k = k + 1;
            end
        end
    end
    if((nom3>=1)||(non3>=1))
        for i = 1:nom3
            for j = 1:non3
                a1 = (i-1)*mt + 1; b1 = i*mt;
                a2 = (j-1)*nt + 1; b2 = j*nt;
                Frs(:,:,1,k) = Fr3(a1:b1,a2:b2,1,ii);
                k = k + 1;
            end
        end
    end
end

end