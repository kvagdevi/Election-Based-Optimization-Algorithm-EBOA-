function [Thresholds,mea,standereddeviation,minimum] = firefly(h,Level,NumberOfIter);
% parameters [n N_iteration alpha betamin gamma]
d=Level;%Dimensions
P=10*d;  %%%population
para=[P NumberOfIter 0.5 0.1 0.1];
% Simple bounds/limits for d-dimensional problems
Lb=1.*ones(1,d);
Ub=255.*ones(1,d);
% for pp=1:1
% Initial random guess
u0=Lb+(Ub-Lb).*rand(1,d);
for pp=1:NumberOfIter
[Thresholds,fval]=ffa_mincon(@cost,u0,Lb,Ub,para,d,h);

% Display results
% bestsolution=u;
 bestojb(pp)=max(fval);
 end;
minimum=max(bestojb)
mea=mean(bestojb)
standereddeviation=std(bestojb)
%%% Put your own cost/objective function here --------%%%
%% Cost or Objective function

% Start FA
function [nbest,fbest]...
           =ffa_mincon(fhandle,u0, Lb, Ub, para,d,h)
% Check input parameters (otherwise set as default values)
if nargin<5, para=[50 500 0.25 0.20 1]; end
if nargin<4, Ub=[]; end
if nargin<3, Lb=[]; end
if nargin<2,
disp('Usuage: FA_mincon(@cost,u0,Lb,Ub,para,cbook)');
end

% n=number of fireflies
% MaxGeneration=number of pseudo time steps
% ------------------------------------------------
% alpha=0.25;      % Randomness 0--1 (highly random)
% betamn=0.20;     % minimum value of beta
% gamma=1;         % Absorption coefficient
% ------------------------------------------------
n=para(1);  MaxGeneration=para(2);
alpha=para(3); betamin=para(4); gamma=para(5);
if length(Lb) ~=length(Ub),
    disp('Simple bounds/limits are improper!');
    return
end

% Calcualte dimension
d=length(u0);

% Initial values of an array
zn=ones(n,1)*10^100;
% ------------------------------------------------
% generating the initial locations of n fireflies
[ns,Lightn]=init_ffa(n,d,Lb,Ub,u0);

% Iterations or pseudo time marching
for k=1:MaxGeneration,     %%%%% start iterations

% This line of reducing alpha is optional
%  alpha=alpha_new(alpha,MaxGeneration);

% Evaluate new solutions (for all n fireflies)
for i=1:n,
   [zn(i)]=fhandle(ns(i,:),h);
   Lightn(i)=zn(i);
end;
   

% Ranking fireflies by their light intensity/objectives
[Lightn,Indexx]=sort(zn,'descend');
ns_tmp=ns;
for i=1:n,
 ns(i,:)=ns_tmp(Indexx(i),:);
end

%% Find the current best
nso=ns; Lighto=Lightn;
nbest=ns(1,:); Lightbest=Lightn(1);

% For output only
fbest(k)=Lightbest;

% Move all fireflies to the better locations
[ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,nbest,...
      Lightbest,alpha,betamin,gamma,Lb,Ub);

end   %%%%% end of iterations
% figure;
% plot(1:1000,fbest);
function [ns,Lightn]=init_ffa(n,d,Lb,Ub,u0)
  % if there are bounds/limits,
if length(Lb)>0,
   for i=1:n,
   ns(i,:)=Lb+(Ub-Lb).*rand(1,d);
   end
else
   for i=1:n,
   ns(i,:)=u0+randn(1,d);
   end
end

Lightn=ones(n,1)*10^100;

% Move all fireflies toward brighter ones
function [ns]=ffa_move(n,d,ns,Lightn,nso,Lighto,...
             nbest,Lightbest,alpha,betamin,gamma,Lb,Ub)
% Scaling of the system
scale=abs(Ub-Lb);
% scale=1;
% Updating fireflies
for i=1:n,
   for j=1:n,
      r=sqrt(sum((ns(i,:)-ns(j,:)).^2));
      % Update moves
if Lightn(i)< Lighto(j), % Brighter and more attractive
   beta0=1; beta=(beta0-betamin)*exp(-gamma*r.^2)+betamin;
   tmpf=alpha.*(rand(1,d)-0.5).*scale;
   ns(i,:)=ns(i,:).*(1-beta)+ns(j,:).*beta+tmpf;
      end
   end % end for j

end % end for i

[ns]=findlimits(n,ns,Lb,Ub);

function [ns]=findlimits(n,ns,Lb,Ub)
for i=1:n,
     % Apply the lower bound
  ns_tmp=ns(i,:);
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);

  % Apply the upper bounds
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move
  ns(i,:)=ns_tmp;
end

%%%%%%%%%%%%%%%%%  Powell Singular 2 Function  -4 to 5 %%%%%%%%
function level = cost(total,histogramCounts)
%% OTSU automatic thresholding method
sumB = 0;
histogramCounts=histogramCounts';
wB = 0;
maximum = 0.0;
threshold1 = 0.0;
threshold2 = 0.0;
sum1 = sum((1:256).*histogramCounts); % the above code is replace with this single line
for ii=1:256
    wB = wB + histogramCounts(ii);
    if (wB == 0)
        continue;
    end
    wF = total - wB;
    if (wF == 0)
        break;
    end
    sumB = sumB +  ii * histogramCounts(ii);
    mB = sumB / wB;
    mF = (sum1 - sumB) ./ wF;
    between = wB .* wF .* (mB - mF) .* (mB - mF);
    if ( between >= maximum )
        threshold1 = ii;
        if ( between > maximum )
            threshold2 = ii;
        end
        maximum = between;
    end
end
level = (threshold1 + threshold2 )/(2);
