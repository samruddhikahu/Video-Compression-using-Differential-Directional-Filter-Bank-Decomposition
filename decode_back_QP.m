function [CurrFrLr, CurrFrar, CurrFrbr] = decode_back_QP(Sband1, Sband2, Sband3, PrevFrL, PrevFra, PrevFrb, MV1, MV2, MV3, QPp)
% Using the decoded frames for ME & MC:-
dfilter = '9-7';

% Re-quantization:-
Sband1 = Sband1.*QPp;
Sband2 = Sband2.*QPp;
Sband3 = Sband3.*QPp;

coeffsr1 = rearrange_back(round(Sband1),3);
coeffsr2 = rearrange_back(round(Sband2),3);
coeffsr3 = rearrange_back(round(Sband3),3);

FDiffr1 = dfbrec(coeffsr1,dfilter);
FDiffr2 = dfbrec(coeffsr2,dfilter);
FDiffr3 = dfbrec(coeffsr3,dfilter);

FComp1r = motionComp(PrevFrL,MV1,16);
FComp2r = motionComp(PrevFra,MV2,8);
FComp3r = motionComp(PrevFrb,MV3,8);

CurrFrLr = FComp1r + FDiffr1;
CurrFrar = FComp2r + FDiffr2;
CurrFrbr = FComp3r + FDiffr3;
end