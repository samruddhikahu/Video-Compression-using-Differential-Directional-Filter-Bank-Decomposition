% This is the code for the compression and decompression of a video.
% Compression is done in CIE L*a*b* color space. It uses EPZS for motion
% estimation and SPECK for entropy coding the Directional Filter Bank
% transformed coefficients.
% The functions dfbdec and dfbrec for Directional Filter Bank decomposition
% are taken from the Contourlet Toolbox by Do and Vetterli.
% Input:- video_mov. Video in .yuv format stored at the same location as this code.
%         The width, height and the chrominance subsampling ratio should
%         also be known beforehand.
% Output:- RecFr. Reconstructed video as a 4D array. MATLAB's built-in
%          implay function is used to play it.
%          MSE, PSNR and SSIM is calculated in CIE L*a*b* space only and
%          stored in a results.txt file along with calculated CR and bpp.
%          Compressed bitstream is represented by arrays bitstream_frs,
%          bitstream_S and also by motion vector difference encoded bits
%          MVD1bits, MVD2bits and MVD3bits. It is stored in the
%          bitstream.txt file.

clc;
clear all;
close all;
Tstrt = tic;

%Obj = VideoReader('vipmosaicking.avi');
%load('tennis_mov.mat');
video_mov = yuv2mov('container_qcif.yuv',176,144,'420');
nframes = length(video_mov);
if(mod(nframes,10)>0)
    nframes = (floor(nframes/10))*10;
end
for i = 1:nframes
    F(:,:,:,i) = video_mov(i).cdata;
end
First = F(:,:,:,1);
[H,W] = size(First(:,:,1,1));
sz = [H W];
TileSize = 8;
save('sz.mat');
%nframes = Obj.NumberOfFrames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENCODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tenc = tic;
tic;
%First = read(Obj,1);
sizes_L = sizes(sz./TileSize,1);
dfilter = '9-7';
% sizes_a = sizes(sz./4,1);
% sizes_b = sizes(sz./4,1);

[h, g] = pfilters1('9-7');
lvl = 3; % Decomposition Level:-
GOPlen = 10;           % GOP Length = 10;
fn = ceil(nframes/GOPlen);    % Number of GOPs 
QPi = 1; QPp = 1;
%k = 1;
mbSize = 16;

% Array initializations:-
MV1 = zeros(2, H*W/mbSize^2, nframes);
MV2 = zeros(2, H*W/mbSize^2, nframes);
MV3 = zeros(2, H*W/mbSize^2, nframes);

MVD1 = zeros(2, H*W/mbSize^2, nframes);
MVD2 = zeros(2, H*W/mbSize^2, nframes);
MVD3 = zeros(2, H*W/mbSize^2, nframes);  

FirstFrL_dwt = zeros(H, W, 1, fn);
FirstFra_dwt = zeros(H/2, W/2, 1, fn);
FirstFrb_dwt = zeros(H/2, W/2, 1, fn);

Sband1 = zeros(H, W, 1, nframes);
Sband2 = zeros(H/2, W/2, 1, nframes);
Sband3 = zeros(H/2, W/2, 1, nframes);

for i = 1:fn
    ff = GOPlen*(i - 1) + 1; %(10*i)-9;
    %First = rgbtoxyzLab(F(:,:,:,ff));
    First = rgb2lab(F(:,:,:,ff));
    FirstFrL = First(:,:,1);
    FirstFra = chr_subsamp(First(:,:,2));
    FirstFrb = chr_subsamp(First(:,:,3));
    
    FirstFrL_dwt(:,:,1,i) = (wfb2dec1(FirstFrL,lvl,h,g))./QPi;
    FirstFra_dwt(:,:,1,i) = (wfb2dec1(FirstFra,lvl,h,g))./QPi;
    FirstFrb_dwt(:,:,1,i) = (wfb2dec1(FirstFrb,lvl,h,g))./QPi;
    
    PrevFrL = wfb2rec1(round(FirstFrL_dwt(:,:,1,i).*QPi),lvl,h,g);%FirstFrL;
    PrevFra = wfb2rec1(round(FirstFra_dwt(:,:,1,i).*QPi),lvl,h,g);%FirstFra; 
    PrevFrb = wfb2rec1(round(FirstFrb_dwt(:,:,1,i).*QPi),lvl,h,g);%FirstFrb; 
    
    %p = 1;
    %MV1f = []; MV2f = []; MV3f = [];
    PrevFrameSAD1 = zeros(1,H*W/mbSize^2);
    PrevFrameSAD2 = zeros(1,H*W/mbSize^2);
    PrevFrameSAD3 = zeros(1,H*W/mbSize^2);
    
    for x = (ff+1):(GOPlen*i)
        % Read the current frame:-
        CurrFr = rgb2lab(F(:,:,:,x));
        %CurrFr = rgbtoxyzLab(F(:,:,:,x));
        CurrFrL = CurrFr(:,:,1);
        CurrFra = chr_subsamp(CurrFr(:,:,2));
        CurrFrb = chr_subsamp(CurrFr(:,:,3));
        
        % Actual difference between the 3 planes of current and previous
        % frames:-
        if(x==2)
            MV1t2 = zeros(2, H*W/mbSize^2);
            MV2t2 = zeros(2, H*W/mbSize^2);
            MV3t2 = zeros(2, H*W/mbSize^2);
        else
            MV1t2 = MV1(:,:,x-2);
            MV2t2 = MV2(:,:,x-2);
            MV3t2 = MV3(:,:,x-2);
        end
        [FComp1, MinFrameSAD1, MV1(:,:,x), MVD1(:,:,x)] = MotionEstComp2(PrevFrL, CurrFrL, PrevFrameSAD1, MV1(:,:,x-1), MV1t2, mbSize, 'EPZS');
        [FComp2, MinFrameSAD2, MV2(:,:,x), MVD2(:,:,x)] = MotionEstComp2(PrevFra, CurrFra, PrevFrameSAD2, MV2(:,:,x-1), MV2t2, mbSize/2, 'EPZS');
        [FComp3, MinFrameSAD3, MV3(:,:,x), MVD3(:,:,x)] = MotionEstComp2(PrevFrb, CurrFrb, PrevFrameSAD3, MV3(:,:,x-1), MV3t2, mbSize/2, 'EPZS');
        
        PrevFrameSAD1 = MinFrameSAD1;
        PrevFrameSAD2 = MinFrameSAD2;
        PrevFrameSAD3 = MinFrameSAD3;
        
        FDiff1 = CurrFrL - FComp1;
        FDiff2 = CurrFra - FComp2;
        FDiff3 = CurrFrb - FComp3;
        
        % Applying Contourlet decomposition on the 3 difference frame planes:-
        coeffs1 = dfbdec(FDiff1,dfilter,3);
        coeffs2 = dfbdec(FDiff2,dfilter,3);
        coeffs3 = dfbdec(FDiff3,dfilter,3);
        
        Sband1(:,:,1,x) = rearrange1(coeffs1,3);
        Sband2(:,:,1,x) = rearrange1(coeffs2,3);
        Sband3(:,:,1,x) = rearrange1(coeffs3,3);
        
        % Quantization in DFB domain:-
        Sband1(:,:,1,x) = round(Sband1(:,:,1,x)./QPp);
        Sband2(:,:,1,x) = round(Sband2(:,:,1,x)./QPp);
        Sband3(:,:,1,x) = round(Sband3(:,:,1,x)./QPp);
        
        % Uisng the decoded frame for ME & MC:-
        [CurrFrLr, CurrFrar, CurrFrbr] = decode_back_QP(Sband1(:,:,1,x),Sband2(:,:,1,x),Sband3(:,:,1,x),PrevFrL,PrevFra,PrevFrb,MV1(:,:,x),MV2(:,:,x),MV3(:,:,x), QPp);
        
        %k = k + 1;
        %p = p + 1;
        
        PrevFrL = CurrFrLr;
        PrevFra = CurrFrar;
        PrevFrb = CurrFrbr;
    end
%     MV1f = cat(3, MV1f, MV1);
%     MV2f = cat(3, MV2f, MV2);
%     MV3f = cat(3, MV3f, MV3);
end
toc

% Encoding first frames using SPECK encoder:-
tic;

sizes_frL = sizes(size(FirstFrL(:,:,1,1))./TileSize,3);
% sizes_fra = sizes(size(FirstFra(:,:,1,1))./2,3);
% sizes_frb = sizes(size(FirstFrb(:,:,1,1))./2,3);

% FirstFrL = round(FirstFrL);
% FirstFra = round(FirstFra);
% FirstFrb = round(FirstFrb);
%H1 = H/2; W1 = W/2;
fid = fopen('bitstream.txt','w+');
fprintf(fid, '%s\r\n', 'bitstream_frs');

FirstFrs = tiled_frames(FirstFrL_dwt,FirstFra_dwt,FirstFrb_dwt,sz./TileSize);

frlen = size(FirstFrs,4);
total_bits_fr = 0;
for i = 1:frlen    
    
    [bp_max_frs(i),bitstream_frs{i},Stat_frs{i}] = SPECK_encoder1_fast2(round(FirstFrs(:,:,1,i)),sizes_frL);    
    total_bits_fr = total_bits_fr + length(bitstream_frs{i});
    fprintf(fid, '%s\r\n', bitstream_frs{i});
end
toc

% Encoding the contourlet sub-bands using SPECK encoder:-
tic;
total_bits = 0;

% Removing the zero frames:-
%[Sband1, Sband2, Sband3] = rem_zero_frames(Sband1, Sband2, Sband3, GOPlen);
Sband1(:,:,:,1:GOPlen:nframes) = [];
Sband2(:,:,:,1:GOPlen:nframes) = [];
Sband3(:,:,:,1:GOPlen:nframes) = [];

Sbands = tiled_frames(Sband1,Sband2,Sband3,sz./TileSize);
Slen = size(Sbands,4);

for i = 1:Slen
    [bp_max_S(i),bitstream_S{i},Stat_S{i}] = SPECK_encoder1_fast2(Sbands(:,:,1,i),sizes_L);
    
    total_bits = total_bits + length(bitstream_S{i});
    fprintf(fid, '%s\r\n', bitstream_S{i});
end
toc

% Encoding the Motion Vector Differences using Exponential Golomb Coding:-
tic;
% Removing the zero MV Differences:-
MVD1(:,:,1:GOPlen:nframes) = [];
MVD2(:,:,1:GOPlen:nframes) = [];
MVD3(:,:,1:GOPlen:nframes) = [];

MVD1bits = MVD_enc(MVD1);
MVD2bits = MVD_enc(MVD2);
MVD3bits = MVD_enc(MVD3);
mv_fr = size(MVD1,2);
nfr = size(MVD1,3);
toc
Tcomp = toc(Tenc)

% Saving the compressed bitstream arrays to a text file:-

fprintf(fid,'%s\r\n %s\r\n %s\r\n %s\r\n %s\r\n %s\r\n %s\r\n %s\r\n %s\r\n','MV1bits', MVD1bits{1}, MVD1bits{2}, 'MV2bits', MVD2bits{1}, MVD2bits{2}, 'MV3bits', MVD3bits{1}, MVD3bits{2});
fprintf(fid,'\r\n');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tdec = tic;

% Decoding the Motion Vector Differences using Exponential Golomb Coding;
tic;
MVD1_rec = MVD_dec(MVD1bits,mv_fr,nfr);
MVD2_rec = MVD_dec(MVD2bits,mv_fr,nfr);
MVD3_rec = MVD_dec(MVD3bits,mv_fr,nfr);
toc

% Decoding first frames using SPECK decoder:-
tic;
lvl = 3; % Decomposition Level.
br_fr = max(bp_max_frs);
for i = 1:frlen
    [FirstFrs_rec(:,:,1,i),FFs_pdwt,Stat1_frs{i}] = SPECK_decoder2_fast2(br_fr,bp_max_frs(i),bitstream_frs{i},sizes_frL);
end
[FirstFrL_dwtr,FirstFra_dwtr,FirstFrb_dwtr] = tiled_frames_rec(FirstFrs_rec,sz./TileSize,[H,W],[H/2,W/2],[H/2,W/2]);
toc

%%% Decoding the contourlet sub-bands using SPECK decoder:-
tic;
br_S = max(bp_max_S);
for i = 1:Slen
    %br = 8;
    [Sbands_rec(:,:,1,i), Sbands_p, Stat1_S{i}] = SPECK_decoder2_fast2(br_S,bp_max_S(i),bitstream_S{i}, sizes_L);
    
end
[Sband1_rec,Sband2_rec,Sband3_rec] = tiled_frames_rec(Sbands_rec,sz./TileSize,[H,W],[H/2,W/2],[H/2,W/2]);
toc

tic;
clear coeffs1 coeffs2 coeffs3;

k = 1; n = 0;
for i = 1:fn
    ff = GOPlen*(i-1) + 1;  %(10*i)-9;
    FirstFrL_rec = wfb2rec1(FirstFrL_dwtr(:,:,1,i).*QPi,lvl,h,g);
    FirstFra_rec = wfb2rec1(FirstFra_dwtr(:,:,1,i).*QPi,lvl,h,g);
    FirstFrb_rec = wfb2rec1(FirstFrb_dwtr(:,:,1,i).*QPi,lvl,h,g);
    
    PrevFrL = FirstFrL_rec;
    PrevFra = FirstFra_rec;
    PrevFrb = FirstFrb_rec;
    First(:,:,1) = FirstFrL_rec;
    First(:,:,2) = chr_resamp(PrevFra);
    First(:,:,3) = chr_resamp(PrevFrb);
    
    n = n + 1;
    RecFr1(:,:,:,n) = First;
    %RecFr(:,:,:,n) = Labtoxyzrgb(First);
    RecFr(:,:,:,n) = lab2rgb(First,'OutputType','uint8');
    
    for x = (ff+1):(GOPlen*i)
        
        % Re-quantization:-
        Sband1_rec(:,:,1,k) = Sband1_rec(:,:,1,k).*QPp;
        Sband2_rec(:,:,1,k) = Sband2_rec(:,:,1,k).*QPp;
        Sband3_rec(:,:,1,k) = Sband3_rec(:,:,1,k).*QPp;

        coeffs1 = rearrange_back(Sband1_rec(:,:,1,k),3);
        coeffs2 = rearrange_back(Sband2_rec(:,:,1,k),3);
        coeffs3 = rearrange_back(Sband3_rec(:,:,1,k),3);
        
        FDiffr1 = dfbrec(coeffs1,dfilter);
        FDiffr2 = dfbrec(coeffs2,dfilter);
        FDiffr3 = dfbrec(coeffs3,dfilter);
        
        FComp1r = motionComp2(PrevFrL, MVD1_rec(:,:,k), mbSize);
        FComp2r = motionComp2(PrevFra, MVD2_rec(:,:,k), mbSize/2);
        FComp3r = motionComp2(PrevFrb, MVD3_rec(:,:,k), mbSize/2);
        k = k + 1;
        
        CurrFrL = FComp1r + FDiffr1;
        CurrFra = FComp2r + FDiffr2;
        CurrFrb = FComp3r + FDiffr3;
        
%         CurrFrL = PrevFrL + FDiffr1;
%         CurrFra = PrevFra + FDiffr2;
%         CurrFrb = PrevFrb + FDiffr3;
        
        CurrFr(:,:,1) = CurrFrL;
        CurrFr(:,:,2) = chr_resamp(CurrFra);
        CurrFr(:,:,3) = chr_resamp(CurrFrb);
        
        PrevFrL = CurrFrL;
        PrevFra = CurrFra;
        PrevFrb = CurrFrb;
        
        n = n + 1;
        RecFr1(:,:,:,n) = CurrFr;
        %RecFr(:,:,:,n) = Labtoxyzrgb(CurrFr);
        RecFr(:,:,:,n) = lab2rgb(CurrFr,'OutputType','uint8');
    end
end
Tdecomp = toc(Tdec)

% Calculating bpp, CR, MSE and PSNR of the video:-
tic;
no_bits = total_bits + total_bits_fr + length(MVD1bits{1}) + length(MVD1bits{2})...
          + length(MVD2bits{1}) + length(MVD2bits{2})...
          + length(MVD3bits{1}) + length(MVD3bits{2});
bpp = no_bits/(H*W*nframes*3);
CR = 8/bpp;

clear Er;
for i = 1:nframes
    clear Vf1 Vf2 Er;
    Vf1 = F(:,:,:,i);
    Vf2 = RecFr(:,:,:,i);
    Er = (Vf1 - Vf2).^2;
    mse(i) = sum(sum(sum(Er)))/(W*H*3);
    psnr(i) = 10*log10(65025/mse(i));
    Ssim(i) = 0;
  for k = 1:3
  Ssim(i) = Ssim(i) + ssim(Vf1(:,:,k), Vf2(:,:,k));
  end
  Ssim(i) = Ssim(i)/3;
end

ovssim = sum(Ssim)/nframes;
ovmse = sum(mse)/nframes;
ovpsnr = 10*log10(65025/ovmse);
toc

% Calculating MSE, PSNR of L, a* and b* planes individually:-
tic;
Ssim1 = 0; Ssim2 = 0; Ssim3 = 0;
for i = 1:nframes
    clear Vf1 Vf2 Er;
    Vf1 = rgb2lab(F(:,:,:,i));
    %Vf1 = rgbtoxyzLab(F(:,:,:,i));
    Vf2 = RecFr1(:,:,:,i);
    Er = (Vf1 - Vf2).^2;
    mse1(i) = sum(sum(Er(:,:,1)))/(W*H);
    mse2(i) = sum(sum(Er(:,:,2)))/(W*H);
    mse3(i) = sum(sum(Er(:,:,3)))/(W*H);
    msex(i) = sum(sum(sum(Er)))/(W*H*3);
    %Ssim(i) = 0;
  %for k = 1:3
    Ssim1 = Ssim1 + ssim(Vf1(:,:,1), Vf2(:,:,1));
    Ssim2 = Ssim2 + ssim(Vf1(:,:,1), Vf2(:,:,1));
    Ssim3 = Ssim3 + ssim(Vf1(:,:,1), Vf2(:,:,1));
  %end
    %Ssim = (Ssim1 + Ssim2 + Ssim3)/3;
  %Ssim_y(i) = ssim(Vf1(:,:,k),Vf2(:,:,k));
end
%ovssimx = sum(Ssim)/nframes;
ovmsex = sum(msex)/nframes;
ovpsnrx = 10*log10(65025/ovmsex);

%ovssim_y = sum(Ssim_y)/nframes;
ovSsim1 = sum(Ssim1)/nframes;
ovSsim2 = sum(Ssim2)/nframes;
ovSsim3 = sum(Ssim3)/nframes;
ovSsimx = (Ssim1 + Ssim2 + Ssim3)/3;

ovmse1 = sum(mse1)/nframes;
ovmse2 = sum(mse2)/nframes;
ovmse3 = sum(mse3)/nframes;

ovpsnr1 = 10*log10(65025/ovmse1);
ovpsnr2 = 10*log10(65025/ovmse2);
ovpsnr3 = 10*log10(65025/ovmse3);
toc

% Saving the CR, MSE, PSNR and SSIM:-
fid = fopen('results.txt','w+');
fprintf(fid,'%12s %12s %12s %12s %12s\r\n','Parameter', 'L*', 'a*', 'b*', 'Overall');
fprintf(fid,'%12s %12.2f %12.2f %12.2f %12.2f\r\n','MSE', ovmse1, ovmse2, ovmse3, ovmsex);
fprintf(fid,'%12s %12.2f %12.2f %12.2f %12.2f\r\n','PSNR', ovpsnr1, ovpsnr2, ovpsnr3, ovpsnrx);
fprintf(fid,'%12s %12.4f %12.4f %12.4f %12.4f\r\n','SSIM', ovSsim1, ovSsim2, ovSsim3, ovSsimx);
fprintf(fid,'\r\n');
fclose(fid);
