function [Fr1,Fr2,Fr3] = tiled_frames_rec(Frs,TileSize,sz1,sz2,sz3)
mt = round(TileSize(1,1)); nt = round(TileSize(1,2));

% m1 = sz1(1,1); n1 = sz1(1,2);
% m2 = sz2(1,1); n2 = sz2(1,2);
% m3 = sz3(1,1); n3 = sz3(1,2);
% if(mt~=TileSize(1,1))
%     x1 = m1/TileSize(1,1);
%     
% end
% if(nt~=TileSize(1,2))
% end

nom1 = sz1(1,1)/TileSize(1,1); non1 = sz1(1,2)/TileSize(1,2);
nom2 = sz2(1,1)/TileSize(1,1); non2 = sz2(1,2)/TileSize(1,2);
nom3 = sz3(1,1)/TileSize(1,1); non3 = sz3(1,2)/TileSize(1,2);

p = 1; ii = 1;
while(ii<=size(Frs,4))
    for i = 1:nom1
        for j = 1:non1
            a1 = (i-1)*mt + 1; b1 = i*mt;
            a2 = (j-1)*nt + 1; b2 = j*nt;
            Fr1(a1:b1,a2:b2,1,p) = Frs(:,:,1,ii);
            ii = ii + 1;
        end
    end
    if((nom2>=1)||(non2>=1))
        for i = 1:nom2
            for j = 1:non2
                a1 = (i-1)*mt + 1; b1 = i*mt;
                a2 = (j-1)*nt + 1; b2 = j*nt;
                Fr2(a1:b1,a2:b2,1,p) = Frs(:,:,1,ii);
                ii = ii + 1;
            end
        end
    end
    if((nom3>=1)||(non3>=1))
        for i = 1:nom3
            for j = 1:non3
                a1 = (i-1)*mt + 1; b1 = i*mt;
                a2 = (j-1)*nt + 1; b2 = j*nt;
                Fr3(a1:b1,a2:b2,1,p) = Frs(:,:,1,ii);
                ii = ii + 1;
            end
        end
    end
    p = p + 1;
end

if(mt~=TileSize(1,1))
    Fr1(1:sz1(1,1),:,1,:) = 0;
    Fr2(1:sz2(1,1),:,1,:) = 0;
    Fr3(1:sz3(1,1),:,1,:) = 0;
end
if(nt~=TileSize(1,2))
    Fr1(:,1:sz1(1,2),1,:) = 0;
    Fr2(:,1:sz2(1,2),1,:) = 0;
    Fr3(:,1:sz3(1,2),1,:) = 0;
end


end