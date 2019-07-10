This repository contains the code for algorithm proposed in the paper titled, "A low-complexity, sequential video compression scheme using
frame differential directional filter bank decomposition in CIE La*b* color space," published in IEEE Access, vol. 5, pp. 14914-14929, 2017.
If using this repository, please cite the above paper.

################################################
HOW TO USE:-
1. There are two versions of the code, one using EPZS for motion estimation and the other one using ARPS for motion estimation.
2. Main compression and decompression codes are contained in the MATLAB files namely 'DDFB_video_compression_EPZS.m' and 'DDFB_video_compression_ARPS.m'
3. Image to be compressed should be in .yuv format stored in the same folder as the code itself.
4. The variables 'QPi' and 'QPp' are the quantization parameters for the intra and inter frames.
5. Apart from these, the parameters 'tilesize' and 'GOPlen' can be changed to control the compression performance of the algorithm.
6. The compressed bitstream is stored in a 'bitstream.txt' file. Compression parameters such as CR, MSE, PSNR and SSIM are stored in a results.txt file.
7. Reconstructed video is represented by a variable 'RecFr' and can be played using the MATLAB's built-in implay function.