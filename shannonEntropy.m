function entropy=shannonEntropy(x,h)
% x=[1 20 20 20 60 80];
% h=rgb2gray(imread('lena.jpg'));
% h=imhist(h);
% h = h / numel(h);
% h(21,1)=-50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS : x is the threshold vector
%          h is the histrogram of the image
% OUTPUT:  entropy is the Entropy of the image after probabilistic 
%          partition by the threshold vector 'h'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by : Sujoy Paul, Jadavpur University, Kolkata %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=round(x);
entropy=0;
x=[1 (x+1) 256];
x(x>256)=256;
s=size(x);
for kk=1:size(x,2)
    if x(kk)<0
        x(kk)=-x(kk);
    elseif x(kk)==0
        x(kk)=x(kk)+1;
    end;
end;
x=sort(x,'ascend');
for i=1:s(2)-1
    temp = h(x(i):x(i+1));
    tsum= sum(temp);
    ent=0;
    if tsum~=0
    for j=x(i):x(i+1)
        if h(j)~=0
            a=h(j)/tsum;
            ent = ent + a*log(a);
        end
    end
    end
    entropy = entropy + ent;    
end

entropy=-entropy;