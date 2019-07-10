function MVp = motion_vector_pred(mVect, imgsz, mbSize, k)

row = imgsz(1,1); col = imgsz(1,2);
nr = row/mbSize; nc = col/mbSize;
%mbc = length(mVect);
MVp = zeros(2,1);


if(k <= nc)
    MVp = mVect(:,k-1); % MV of A block. MVleft.
elseif(mod(k-1,nc)==0)
    MVp(1,1) = ceil(median([mVect(1, k - nc), mVect(1, k - nc + 1)])); %MV of B and C blocks. MVtop and MVtop-right.
    MVp(2,1) = ceil(median([mVect(2, k - nc), mVect(2, k - nc + 1)]));
elseif(mod(k,nc)==0)
    MVp(1,1) = ceil(median([mVect(1, k - 1), mVect(1, k - nc), mVect(1, k - nc - 1)])); %MV of A, B and D blocks. MVtop and MVtop-left and MVleft.
    MVp(2,1) = ceil(median([mVect(2, k - 1), mVect(2, k - nc), mVect(2, k - nc - 1)]));
else
    MVp(1,1) = ceil(median([mVect(1, k - 1), mVect(1, k - nc), mVect(1, k - nc + 1)])); %MV of A, B and C blocks. MVtop and MVtop-right and MVleft.
    MVp(2,1) = ceil(median([mVect(2, k - 1), mVect(2, k - nc), mVect(2, k - nc + 1)]));
end
end