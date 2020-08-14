function [nvox, QQ, targets] = se_defineVoxel(xyzmm,locNr,mode)
global st;
global MAP;
global SPM;
global xSPM;
global CLUSTER;
global displayType;
global group
global defaults;


if mode == 2
    QQ = (MAP(locNr).LR == sign(xyzmm(1)));
    xyz = group(1).xSPM.xVOL.iM*[MAP(locNr).XYZmm(:,QQ);ones(1,size(QQ(QQ),2))];
    QQ(QQ) = spm_sample_vol(spm_vol([SPM.swd '/mask.img']),xyz(1,:),xyz(2,:),xyz(3,:),0)>0;
else
    dec = 0; QQ = zeros(size(MAP(locNr).Z)); mx = max(  MAP(locNr).Z(  sign(MAP(locNr).XYZmm(1,:)) == sign(xyzmm(1))  )  );
    % def PosVox
    sideMask = (MAP(locNr).LR == sign(xyzmm(1)));
    xyz = group(1).xSPM.xVOL.iM*[MAP(locNr).XYZmm(:,sideMask);ones(1,size(MAP(locNr).XYZmm(:,sideMask),2))];
    posVox = min([prod(size(find(MAP(locNr).Z(:,sideMask))))/10 size(find(spm_sample_vol(spm_vol([SPM.swd '/mask.img']),xyz(1,:),xyz(2,:),xyz(3,:),0)>0),2)]);
    while prod(size(find(QQ))) < posVox
        QQ = (MAP(locNr).Z >= mx-(25*dec)) & sideMask;
        xyz = group(1).xSPM.xVOL.iM*[MAP(locNr).XYZmm(:,QQ);ones(1,size(QQ(QQ),2))];
        QQ(QQ) = spm_sample_vol(spm_vol([SPM.swd '/mask.img']),xyz(1,:),xyz(2,:),xyz(3,:),0)>0;
        dec = dec+1;
        mxx(dec) = prod(size(find(QQ)));
    end           
    dec = dec-1;
	if abs(mxx(dec+1) - posVox) > abs(mxx(dec) - posVox) || dec == 8
        dec = dec-1;
        QQ = (MAP(locNr).Z >= mx-(25*dec)) & sideMask;
        xyz = group(1).xSPM.xVOL.iM*[MAP(locNr).XYZmm(:,QQ);ones(1,size(QQ(QQ),2))];
%        QQ(QQ) = spm_sample_vol(spm_vol([SPM.swd '/mask.img']),xyz(1,:),xyz(2,:),xyz(3,:),0)>0;
	end
end
targets = MAP(locNr).XYZmm(:,QQ);
nvox = prod(size(find(QQ)));