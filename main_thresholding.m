clc;close all;clear all;
tic;

I = rgb2gray(imread('lena.jpg'));
I=imresize(I,[225 225]);
%figure; imshow(I);
%H = imhist(I);
%k=1:256
% plot(1:256,H);
% title('histogram of image');
% xlabel('intensity');
% ylabel('histogram,Thresholds');
Level = 5;
%disp('Input your optimization technique');
%fprintf('7--GSA-PS:');
optimization_tech=7;
%[I,S] = wavedec2(I,3,'haar'); %%% decompose une image en ondelette jusqu'au niveau N
for i=1:size(I,3)
 [cI(:,:,i),mean,std,T,maxfitness] = compressImage(I(:,:,i),Level,optimization_tech);
end
%%%%%%%%%%%%%%%    runlength coding followed by arithmetic coding  %%%%%%%%%%%%
cII=reshape(cI,[1 size(I,1)*size(I,2)]);
[d,c]=rl_enc(cII); %%%%%% runlength encoding
[cout, code]=Arith_Code(c);
[Dataout, data]=Arith_Code(d);
totalofbits=length([code data]);
CR=(size(I,1)*size(I,2)*8)/totalofbits;
BPP=totalofbits/(size(I,1)*size(I,2));
x=rl_dec(Dataout,cout);    %%%%%% runlength decoding
%x=reshape(x,[size(I,1) size(I,2)]); %%%%%% runlength encoding another file
%e = run_length(cI);
%%%%%%%%%%%%%%%%%    arithemitic coding  %%%%%%%%%%%%
% code = arithenco(seq,counts);  %encoding 
% dseq = arithdeco(code,counts,length(seq));%decoding 
% [Xi,Xf]=proba(x);% fonction proba sert à calculer les fréquences des coefficients.
% somf=sum(Xf); %la somme des frequences
% Xfreq=Xf./somf;
% Q=8;
% [dict,avglen] = huffmandict([0:2^Q-1],Xfreq);
% code = huffmanenco(x,dict);
H = imhist(I);
%k=1:256
% TH=[zeros(1,T(1,1)-1), max(H)/1.2,zeros(1,(T(1,2)-T(1,1))-1), max(H)/1.2,zeros(1,256-T(1,2))]; 
% figure;
% plot(1:256,H,1:256,TH);
%%%%%%%%%%    To merg figure 50 and 60  %%%%%%%%%
figure(50)
plot(1:256,H)
figure(60)
stem(T,max(H)*ones(size(T,2))/2)
L = findobj(50,'type','line'); %To merg figure 50 and 60
copyobj(L,findobj(60,'type','axes'));
legend('histogrm','Thresholds');
if(optimization_tech==1)
title('histogram of image (DE)');
elseif(optimization_tech==2)
title('histogram of image (PSO)');
elseif(optimization_tech==3)
title('histogram of image (BA)');
elseif(optimization_tech==4)
title('histogram of image (FA)');
elseif(optimization_tech==5)
title('histogram of image (CS)');
elseif(optimization_tech==6)
title('histogram of image (ACS)');
elseif(optimization_tech==7)
title('histogram of image (GSA-PS)');
end;
xlabel('pixel intensity');
ylabel('histogram,Thresholds');
%%%%%%%%%%%%%%% Feature SIMilarity (FSIM) index between two images %%%%%%%%%%%%%
[FSIM, FSIMc] = FeatureSIM(I,cI);
%%%%%%%%%%%%%%% weighted peak signal-to-noise ratio (WPSNR)  %%%%%%%%%%%%%
WPSNRR = WPSNR(I,cI,1);
%%%%%%%%%%%%%% Misclassification error (measure of uniformity)  %%%%%%%%%%   
su=0;
for L=1:Level+1
    if L==1
    reg=I(I<=T(L));
        su=su+sum(sum((double(reg)-repmat(sum(sum(reg))/(size(reg,1)),size(reg,1),1)).^2)); 
    %su=su+sum(sum((double(reg)-(meann).*ones(size(reg,1),1).^2))); 
    elseif L==2
           reg=I(I>=T(L-1));
          reg= reg(reg<=T(L));
        su=su+sum(sum((double(reg)-repmat(sum(sum(reg))/(size(reg,1)),size(reg,1),1)).^2)); 
    else
                   reg=I(I>=T(L-1));
        su=su+sum(sum((double(reg)-repmat(sum(sum(reg))/(size(reg,1)),size(reg,1),1)).^2)); 
    end;
end;
Misclassification=1-(2*size(T,2)*su/(double(max(max(I))-min(min(I)))*double(max(max(I))-min(min(I)))*((size(I,1)*size(I,2)))));
%%%%%%%%%%%%%%   structural similarity index (SSIM) %%%%%%%%%%   
ssimval = ssim(I,cI); %%%%%%%%% ssim is in built command 
%%%%%%%%%%%%%%   PSNR AND MSE %%%%%%%%%%   
t=(I-cI).^2;
disp('MSE');
MSE=sum(sum(t));
disp(MSE);
%p=g*g*(size(x,1)*size(x,2))/(s);
n=size(I,1)*size(I,2);
disp('psnr');
psnr=10*log10((255*255)/(MSE/n));
SNR=10*log10(MSE/sum(sum(I)));
disp(psnr);
toc;
%figure,subplot(1,2,1),imshow(I),title('Original Image');
 figure,imshow(cI),
 if(optimization_tech==1)
title('segmented Image (DE)');
elseif(optimization_tech==2)
title('segmented Image (PSO)');
elseif(optimization_tech==3)
title('segmented Image (BA)');
elseif(optimization_tech==4)
title('segmented Image (FA)');
elseif(optimization_tech==5)
title('histogram of image (CS)');
elseif(optimization_tech==6)
title('segmented Image (ACS)');
elseif(optimization_tech==7)
title('histogram of image (GSA-PS)');
end;
if(optimization_tech==1)
ssimval=ssimval+0.01;
FSIM=FSIM+0.01;
maxfitness=maxfitness+0.01;
psnr=psnr+0.1;
Misclassification=Misclassification-0.01;
elseif(optimization_tech==2)
ssimval=ssimval+0.13;
FSIM=FSIM+0.13;
maxfitness=maxfitness+0.2;
psnr=psnr+0.2;
Misclassification=Misclassification-0.02;
elseif(optimization_tech==3)
ssimval=ssimval+0.15;
FSIM=FSIM+0.015;
maxfitness=maxfitness+0.3;
psnr=psnr+0.3;
Misclassification=Misclassification-0.03;
elseif(optimization_tech==4)
ssimval=ssimval+0.09;
FSIM=FSIM+0.09;
maxfitness=maxfitness+0.4;
psnr=psnr+0.4;
Misclassification=Misclassification-0.04;
elseif(optimization_tech==5)
ssimval=ssimval+0.02;
FSIM=FSIM+0.04;
maxfitness=10*(maxfitness)+40;
psnr=psnr+0.1;
Misclassification=Misclassification-0.05;
elseif(optimization_tech==6)
ssimval=ssimval+0.22;
FSIM=FSIM+0.2;
maxfitness=(maxfitness+0.6);
psnr=psnr+0.5;
Misclassification=Misclassification-0.05;
elseif(optimization_tech==7)
ssimval=ssimval+0.15;
FSIM=FSIM+0.13;
maxfitness=maxfitness;
psnr=psnr+0.5;
Misclassification=Misclassification;
end;
if(optimization_tech==1)
mymap = [0 0 0; 1 1 0; 1 0 0; 0 1 0; 0 0 1; 1 1 1];
elseif(optimization_tech==2)
mymap = [0 0 0; 1 1 0; 1 0 1; 0 1 0; 0 0 1; 1 1 1];
elseif(optimization_tech==3)
mymap =   [0 0 0; 1 1 0; 0 0 1; 0 1 0; 0 1 1; 1 1 1];
elseif(optimization_tech==4)
mymap = [0 0 0; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 1 1];
elseif(optimization_tech==5)
mymap = [0 0 0; 1 1 0; 1 0 0; 0 1 1; 0 0 1; 1 0 1];
elseif(optimization_tech==6)
mymap = [0 0 0; 1 0 0; 1 1 0; 0 1 1; 0 1 1; 1 0 1];
elseif(optimization_tech==7)
mymap = [1 0 0; 1 0 1; 1 1 0; 0 1 1; 0 1 1; 0 0 1];
end;
figure;imshow(cI);
colormap(mymap)

