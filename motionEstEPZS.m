function [motionVect, MinFrameSAD, EPZScomputations] = motionEstEPZS(imgP,imgI,PrevSAD, mVect1,mVect2,mbSize,p)
% A function to calculate motion vectors using EPZS technique.
[row, col] = size(imgI);
vectors = zeros(2,row*col/mbSize^2);
MinSAD = zeros(1,row*col/mbSize^2);
computations = 0;
mbCount = 1;

for i = 1:mbSize:row-mbSize+1
    for j = 1:mbSize:col-mbSize+1
        checkMatrix = zeros(2*p+1,2*p+1);
        costs = ones(2*p+1,2*p+1)*65537;
        
        % Subset A:-
        if((mbCount>1)&&(mbCount<=(col/mbSize)))
            Med_MV = [vectors(1,mbCount-1), vectors(2,mbCount-1)];
        elseif(mbCount==(col/mbSize)+1)
            Med_MV = [ceil(median([vectors(1,mbCount-1),vectors(1,mbCount-(col/mbSize)+1)])), ceil(median([vectors(2,mbCount-1),vectors(2,mbCount-(col/mbSize)+1)]))];
        elseif(mbCount>(col/mbSize)+1)
            Med_MV = [ceil(median([vectors(1,mbCount-1),vectors(1,mbCount-(col/mbSize)),vectors(1,mbCount-(col/mbSize)+1)])),...
                ceil(median([vectors(2,mbCount-1),vectors(2,mbCount-(col/mbSize)),vectors(2,mbCount-(col/mbSize)+1)]))];
        end
        if(exist('Med_MV'))
            % Calculate SAD at Med_MV.
            sx = i+Med_MV(1,1); sy = j + Med_MV(1,2);
            if ( sx < 1 || sx+mbSize-1 > row || sy < 1 || sy+mbSize-1 > col)
            else
                SAD_med = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1),imgI(sx:sx+mbSize-1,sy:sy+mbSize-1),mbSize);
                checkMatrix(Med_MV(1,1) + p+1, Med_MV(1,2) + p+1) = 1;
                costs(Med_MV(1,1) + p+1, Med_MV(1,2) + p+1) = SAD_med;
                % Check with the threshold.
                if(SAD_med < (256/mbSize^2))
                    vectors(1:2,mbCount) = Med_MV;
                    MinSAD(1,mbCount) = SAD_med;
                    mbCount = mbCount + 1;
                    continue;
                end
            end
        end
        %else
        % Subset B:-
        MV_subB(1,:) = [0,0];
        %MV_subB(2,:) = [mVect1(1,mbCount), mVect1(2,mbCount)]; % MV of collocated MB
        if(mbCount==1)
            nmb = 1;
        elseif((mbCount>1)&&(mbCount<=(col/mbSize)))
            MV_subB(2,:) = [vectors(1,mbCount-1), vectors(2,mbCount-1)]; % left MB
            nmb = 2;
        elseif(mbCount==(col/mbSize)+1)
            MV_subB(2,:) = [vectors(1,mbCount-1), vectors(2,mbCount-1)]; % left MB
            MV_subB(3,:) = [vectors(1,mbCount-(col/mbSize)+1), vectors(2,mbCount-(col/mbSize)+1)]; % top-right MB
            nmb = 3;
        else
            MV_subB(2,:) = [vectors(1,mbCount-1), vectors(2,mbCount-1)]; % left MB
            MV_subB(3,:) = [vectors(1,mbCount-(col/mbSize)+1), vectors(2,mbCount-(col/mbSize)+1)]; % top-right MB
            MV_subB(4,:) = [vectors(1,mbCount-(col/mbSize)), vectors(2,mbCount-(col/mbSize))]; % top MB
            nmb = 4;
        end
        %             MV_subB(1,:) = [0,0];
        %             MV_subB(2,:) = [vectors(1,mbCount-1), vectors(2,mbCount-1)]; % left MB
        %             MV_subB(3,:) = [vectors(1,mbCount-(col/mbSize)), vectors(2,mbCount-(col/mbSize))]; % top MB
        %             MV_subB(4,:) = [vectors(1,mbCount-(col/mbSize)+1), vectors(2,mbCount-(col/mbSize)+1)]; % top-right MB
        %             MV_subB(5,:) = [mVect1(1,mbCount), mVect1(2,mbCount)]; % collocated MB
        %costs = ones(1, 5) * 65537;
        
        for k = 1:nmb
            refBlkVer = i + MV_subB(k,1);   % row/Vert co-ordinate for ref block
            refBlkHor = j + MV_subB(k,2);   % col/Horizontal co-ordinate
            if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                continue;
            end
            
            if (refBlkHor < j-p || refBlkHor > j+p || refBlkVer < i-p || refBlkVer > i+p)
                continue;
            elseif (checkMatrix(MV_subB(k,1)+p+1 , MV_subB(k,2)+p+1) == 1)
                compcost(k) = costs(MV_subB(k,1) + p+1, MV_subB(k,2) + p+1);
                %costs(MV_subB(k,1) + p+1, MV_subB(k,2) + p+1) = 1;
                continue
            end
            
            costs(MV_subB(k,1) + p+1, MV_subB(k,2) + p+1) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), imgI(refBlkVer:refBlkVer+mbSize-1,refBlkHor:refBlkHor+mbSize-1), mbSize);
            checkMatrix(MV_subB(k,1)+p+1, MV_subB(k,2)+p+1) = 1;
            compcost(k) = costs(MV_subB(k,1) + p+1, MV_subB(k,2) + p+1);
            computations =  computations + 1;
        end
        
        SADminB = min(min(costs));
        [point(1,1), point(1,2)] = find(costs==SADminB);
            
        % Calculate threshold T2/T3:- edit
        SADm = min(compcost(2:end));
        if(isempty(SADm))
            SADm = 0;
        end
        SADm = min(SADm,PrevSAD(1,mbCount));
        T2 = (1.2*SADm + 128)/(mbSize^2);
            
        if(SADminB<T2)
            vectors(1:2,mbCount) = [point(1,1)-p-1,point(1,2)-p-1]';
            MinSAD(1,mbCount) = SADminB;
            mbCount = mbCount + 1;
            continue;
        end
        
        % Subset C:-
        MV_subC(1,:) = [mVect1(1,mbCount), mVect1(2,mbCount)]; % MV of collocated MB
        MV_subC(2,:) = [2*mVect1(1,mbCount)-mVect2(1,mbCount), 2*mVect1(2,mbCount)-mVect2(2,mbCount)]; % MVaccelpred;
        if((i==1)&&(j==1))
            MV_subC(3,:) = [mVect1(1,mbCount+1), mVect1(2,mbCount+1)]; % MVrightco;
            MV_subC(4,:) = [mVect1(1,mbCount+(col/mbSize)), mVect1(2,mbCount+(col/mbSize))]; % MVdownco;
            nmc = 4;
        elseif((i==1)&&(j>1)&&(j<(col - mbSize + 1)))
            MV_subC(3,:) = [mVect1(1,mbCount+1), mVect1(2,mbCount+1)]; % MVrightco;
            MV_subC(4,:) = [mVect1(1,mbCount+(col/mbSize)), mVect1(2,mbCount+(col/mbSize))]; % MVdownco;
            MV_subC(5,:) = [mVect1(1,mbCount-1), mVect1(2,mbCount-1)]; % MVleftco;
            nmc = 5;
        elseif((i>1)&&(j==1)&&(i<(row - mbSize + 1)))
            MV_subC(3,:) = [mVect1(1,mbCount+1), mVect1(2,mbCount+1)]; % MVrightco;
            MV_subC(4,:) = [mVect1(1,mbCount+(col/mbSize)), mVect1(2,mbCount+(col/mbSize))]; % MVdownco;
            MV_subC(5,:) = [mVect1(1,mbCount-(col/mbSize)), mVect1(2,mbCount-(col/mbSize))]; % MVtopco;
            nmc = 5;
        elseif((j==1)&&(i==(row - mbSize + 1)))
            MV_subC(3,:) = [mVect1(1,mbCount+1), mVect1(2,mbCount+1)]; % MVrightco;
            MV_subC(4,:) = [mVect1(1,mbCount-(col/mbSize)), mVect1(2,mbCount-(col/mbSize))]; % MVtopco;
            nmc = 4;
        elseif((i==1)&&(j==(col - mbSize + 1)))
            MV_subC(3,:) = [mVect1(1,mbCount-1), mVect1(2,mbCount-1)]; % MVleftco;
            MV_subC(4,:) = [mVect1(1,mbCount+(col/mbSize)), mVect1(2,mbCount+(col/mbSize))]; % MVdownco;
            nmc = 4;
        elseif((i==(row - mbSize + 1))&&(j==(col - mbSize + 1)))
            MV_subC(3,:) = [mVect1(1,mbCount-1), mVect1(2,mbCount-1)]; % MVleftco;
            MV_subC(4,:) = [mVect1(1,mbCount-(col/mbSize)), mVect1(2,mbCount-(col/mbSize))]; % MVtopco;
            nmc = 4;
        elseif((i==(row - mbSize + 1))&&(j>1)&&(j<(col - mbSize + 1)))
            MV_subC(3,:) = [mVect1(1,mbCount-1), mVect1(2,mbCount-1)]; % MVleftco;
            MV_subC(4,:) = [mVect1(1,mbCount+1), mVect1(2,mbCount+1)]; % MVrightco;
            MV_subC(5,:) = [mVect1(1,mbCount-(col/mbSize)), mVect1(2,mbCount-(col/mbSize))]; % MVtopco;
            nmc = 5;
        elseif((j==(col - mbSize + 1))&&(i>1)&&(i<(row - mbSize + 1)))
            MV_subC(3,:) = [mVect1(1,mbCount-1), mVect1(2,mbCount-1)]; % MVleftco;
            MV_subC(4,:) = [mVect1(1,mbCount-(col/mbSize)), mVect1(2,mbCount-(col/mbSize))]; % MVtopco;
            MV_subC(5,:) = [mVect1(1,mbCount+(col/mbSize)), mVect1(2,mbCount+(col/mbSize))]; % MVdownco;
            nmc = 5;
        else
            MV_subC(3,:) = [mVect1(1,mbCount-1), mVect1(2,mbCount-1)]; % MVleftco;
            MV_subC(4,:) = [mVect1(1,mbCount+1), mVect1(2,mbCount+1)]; % MVrightco;
            MV_subC(5,:) = [mVect1(1,mbCount-(col/mbSize)), mVect1(2,mbCount-(col/mbSize))]; % MVtopco;
            MV_subC(6,:) = [mVect1(1,mbCount+(col/mbSize)), mVect1(2,mbCount+(col/mbSize))]; % MVdownco;
            nmc = 6;
        end
        %costs = ones(1, 5) * 65537;
        
        for k = 1:nmc
            refBlkVer = i + MV_subC(k,1);   % row/Vert co-ordinate for ref block
            refBlkHor = j + MV_subC(k,2);   % col/Horizontal co-ordinate
            if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                continue;
            end
            
            if (refBlkHor < j-p || refBlkHor > j+p || refBlkVer < i-p || refBlkVer > i+p)
                continue;
            elseif (checkMatrix(MV_subC(k,1)+p+1 , MV_subC(k,2)+p+1) == 1)
                continue
            end
            
            costs(MV_subC(k,1)+p+1 , MV_subC(k,2)+p+1) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), imgI(refBlkVer:refBlkVer+mbSize-1,refBlkHor:refBlkHor+mbSize-1), mbSize);
            checkMatrix(MV_subC(k,1)+p+1, MV_subC(k,2)+p+1) = 1;
            computations =  computations + 1;
        end
        
        SADminC = min(min(costs));
        [point(1,1), point(1,2)] = find(costs==SADminC);
                
        if(SADminC<T2)
            vectors(1:2,mbCount) = [point(1,1)-p-1,point(1,2)-p-1]';
            MinSAD(1,mbCount) = SADminC;
            mbCount = mbCount + 1;
            continue;
        end
        
        % Continue by using refinement pattern:-
        x = i + point(1,1)-p-1;
        y = j + point(1,2)-p-1;
        
        %costs(3) = SADminC;
        
        % The index points for Small Diamond Search Pattern
        SDSP(1,:) = [ 0 -1];
        SDSP(2,:) = [-1  0];
        SDSP(3,:) = [ 0  0];
        SDSP(4,:) = [ 1  0];
        SDSP(5,:) = [ 0  1];
                    
        costs = ones(1,5)*65537;
        costs(3) = SADminC;
        doneFlag = 0;
        while (doneFlag == 0)
            for k = 1:5
                refBlkVer = x + SDSP(k,1);   % row/Vert co-ordinate for ref block
                refBlkHor = y + SDSP(k,2);   % col/Horizontal co-ordinate
                if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                    continue;
                end
                
                if (k == 3)
                    continue;
%                 elseif (refBlkHor < j-p || refBlkHor > j+p || refBlkVer < i-p ...
%                         || refBlkVer > i+p)
%                     continue;
                end
                if(((x-i+SDSP(k,1)+p+1)>=(2*p + 1))||((y-j+SDSP(k,2)+p+1)>=(2*p + 1))||((x-i+SDSP(k,1)+p+1)<1)||((y-j+SDSP(k,2)+p+1)<1))
                    doneFlag = 1;
                    break;
                elseif(checkMatrix(x-i+SDSP(k,1)+p+1 , y-j+SDSP(k,2)+p+1) == 1)
                    continue;
                end
                
                %    costs(y-i+SDSP(k,2)+p+1 , x-j+SDSP(k,1)+p+1) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                %    imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                costs(k) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                checkMatrix(x-i+SDSP(k,1)+p+1, y-j+SDSP(k,2)+p+1) = 1;
                computations =  computations + 1;
                
            end
            
            [cost,point] = min(costs);
            %[point(1,1), point(1,2)] = find(costs==cost);
            
            %if ((point(1,1) == p+1)&&(point(1,2) == p+1)||(length(find(costs==65537))==15*15))
            if ((point == 3)||(length(find(costs==65537))==5))
                doneFlag = 1;
            else
                x = x + SDSP(point, 1);
                y = y + SDSP(point, 2);
                costs = ones(1,5)*65537;
                costs(3) = cost;
                %costs = ones(2*p+1,2*p+1) * 65537;
                %costs(p+1,p+1) = cost;
            end
            
        end  % while loop ends here
        vectors(1,mbCount) = x - i;    % row co-ordinate for the vector
        vectors(2,mbCount) = y - j;    % col co-ordinate for the vector
        MinSAD(1,mbCount) = cost;
        mbCount = mbCount + 1;
    end
end
MinFrameSAD = MinSAD;
motionVect = vectors;
EPZScomputations = computations/(mbCount-1) ;
end
