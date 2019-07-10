function [FrameComp, MinFrameSAD, motionVect, MVD] = MotionEstComp2(Prev, Curr, PrevSAD, mVect1, mVect2, mbSize, ME_type)
p = 32;
imgI = Prev;
imgP = Curr;
[m, n] = size(imgI);
if(mod(m,mbSize)>0)
    tmp = mod(m,mbSize);
    imgI = [imgI;
        zeros(mbSize-tmp,n)];
    imgP = [imgP;
        zeros(mbSize-tmp,n)];
end
[m, n] = size(imgI);
if(mod(n,mbSize)>0)
    tmp = mod(n,mbSize);
    imgI = [imgI zeros(m,mbSize-tmp)];
    imgP = [imgP zeros(m,mbSize-tmp)];
end

switch ME_type
    case {'ES'}
        % Exhaustive Search
        [motionVect, computations] = motionEstES(imgP,imgI,mbSize,p);
        FrameComp = motionComp(imgI, motionVect, mbSize);
        % ESpsnr = imgPSNR(imgP, FrameComp, 255);
        % EScomputations = computations;
        t = computations;
    case {'TSS'}
        % Three Step Search
        [motionVect,computations ] = motionEstTSS(imgP,imgI,mbSize,p);
        FrameComp = motionComp(imgI, motionVect, mbSize);
        % TSSpsnr = imgPSNR(imgP, FrameComp, 255);
        % TSScomputations = computations;
        t = computations;
    case {'SETSS'}
        % Simple and Efficient Three Step Search
        [motionVect, computations] = motionEstSESTSS(imgP,imgI,mbSize,p);
        FrameComp = motionComp(imgI, motionVect, mbSize);
        % SESTSSpsnr = imgPSNR(imgP, FrameComp, 255);
        % SESTSScomputations = computations;
        t = computations;
    case {'NTSS'}
        % New Three Step Search
        [motionVect,computations ] = motionEstNTSS(imgP,imgI,mbSize,p);
        FrameComp = motionComp(imgI, motionVect, mbSize);
        % NTSSpsnr(i+1) = imgPSNR(imgP, FrameComp, 255);
        % NTSScomputations(i+1) = computations;
        t = computations;
    case {'FSS'}
        % Four Step Search
        [motionVect, computations] = motionEst4SS(imgP,imgI,mbSize,p);
        FrameComp = motionComp(imgI, motionVect, mbSize);
        % SS4psnr = imgPSNR(imgP, FrameComp, 255);
        % SS4computations = computations;
        t = computations;
    case {'DS'}
        % Diamond Search
        [motionVect, computations] = motionEstDS(imgP,imgI,mbSize,p);
        FrameComp = motionComp(imgI, motionVect, mbSize);
        % DSpsnr(i+1) = imgPSNR(imgP, imgComp, 255);
        % DScomputations(i+1) = computations;
        t = computations;
    case {'ARPS'}
        % Adaptive Rood Patern Search
        [motionVect, computations] = motionEstARPS(imgP,imgI,mbSize,p);
        FrameComp = motionComp(imgI, motionVect, mbSize);
        % ARPSpsnr(i+1) = imgPSNR(imgP, imgComp, 255);
        % ARPScomputations(i+1) = computations;
        t = computations;
    case {'EPZS'}
        % Enhanced Predictive Zonal Search
        %mVect1 = MVprev1; mVect2 = MVprev2;
        [motionVect, MinFrameSAD, computations] = motionEstEPZS(imgP,imgI,PrevSAD, mVect1,mVect2,mbSize,p);
        FrameComp = motionComp(imgI, motionVect, mbSize);
        % EPZSpsnr(i+1) = imgPSNR(imgP, imgComp, 255);
        % EPZScomputations(i+1) = computations;
        t = computations;
end

MVD = motion_vector_difference(motionVect, [m, n], mbSize);
%MVD = motionVect - MVp;
end