function sizex = sizes(sz,l)
sizex = ceil(sz);
for i = 1:l
    clear sz1;
    sz1 = ceil(sz./(2^i));
    sizex = [sz1; sizex];
end
%sizex = [sizex(1,:); sizex];
end