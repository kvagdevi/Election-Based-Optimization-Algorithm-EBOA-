function [J, meann, stdd,Thresholds,maxfitness]=compressImage(I,Level,optimization_tech)

if length(size(I))>2
    I=rgb2gray(I);
end
H = imhist(I);
h = H / numel(I);
NumberOfIter = 30;
if optimization_tech==1
    [Thresholds,meann,stdd,maxfitness] = diffEvolution(h,Level,NumberOfIter);
elseif optimization_tech==2
    [Thresholds,meann,stdd,maxfitness] = pso(h,Level,NumberOfIter);
elseif optimization_tech==3
     [Thresholds,meann,stdd,maxfitness] = bat(h,Level,NumberOfIter);
elseif optimization_tech==4
    [Thresholds,meann,stdd,maxfitness] = firefly(h,Level,NumberOfIter);
    elseif optimization_tech==5
    [Thresholds,meann,stdd,maxfitness] = CS(h,Level,NumberOfIter);
elseif optimization_tech==6
    [Thresholds,meann,stdd,maxfitness] = ACS(h,Level,NumberOfIter);
    elseif optimization_tech==7
    [Thresholds,meann,stdd,maxfitness] = GSAPS(h,Level,NumberOfIter);
end;
Thresholds=sort(Thresholds,'ascend');
Thres=[1 Thresholds+1 256];
J=0;
gray_level = zeros(1,length(Thres)-1);

for l=2:length(Thres)
    if sum(h(Thres(l-1):Thres(l)))~=0;
        gray_level(l-1)=sum((h(Thres(l-1):Thres(l)))'.*double(Thres(l-1):Thres(l)))/sum(h(Thres(l-1):Thres(l)));
        J=J+gray_level(l-1)*(I>Thres(l-1) & I<=Thres(l));
    end
end

J=uint8(round(J));
% figure;
% imshow(J)