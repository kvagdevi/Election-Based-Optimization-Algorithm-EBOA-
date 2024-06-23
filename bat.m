function [best,mea,standereddeviation,maxfitness] = bat(h,Level,NumberOfIter);
% x=[1 20 20 20 60 80];
% h=rgb2gray(imread('lena.jpg'));
% h=imhist(h);
% h = h / numel(h);
% h(21,1)=-50;
para=[10*Level NumberOfIter 0.5 1];  
n=para(1);      % Population size, typically 10 to 40
N_gen=para(2);  % Number of generations
A=para(3);      % Loudness  (constant or decreasing)
r=para(4);      % Pulse rate (constant or decreasing)
% This frequency range determines the scalings
% You should change these values if necessary
Qmin=0;         % Frequency minimum
Qmax=30;         % Frequency maximum
% Iteration parameters
N_iter=0;       % Total number of function evaluations
% Dimension of the search variables
d=Level;           % Number of dimensions 
% Lower limit/bounds/ a vector
Lb=0*ones(1,d);
% Upper limit/bounds/ a vector
Ub=255*ones(1,d);
% Initializing arrays
Q=zeros(n,1);   % Frequency
v=zeros(n,d);   % Velocities
% Initialize the population/solutions
for pp=1:NumberOfIter
for i=1:n,
  Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);
  Fitness(i)=shannonEntropy(Sol(i,:),h);
end
% Find the initial best solution
[fmin,I]=max(Fitness);
best=Sol(I,:);


% ======================================================  %
% Note: As this is a demo, here we did not implement the  %
% reduction of loudness and increase of emission rates.   %
% Interested readers can do some parametric studies       %
% and also implementation various changes of A and r etc  %
% ======================================================  %

% Start the iterations -- Bat Algorithm (essential part)  %
for t=1:N_gen, 
% Loop over all bats/solutions
        for i=1:n,
          Q(i)=Qmin+(Qmin-Qmax)*rand;
          v(i,:)=v(i,:)+(Sol(i,:)-best)*Q(i);
          S(i,:)=Sol(i,:)+v(i,:);
          % Apply simple bounds/limits
          Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
          % Pulse rate
          if rand>r
          % The factor 0.001 limits the step sizes of random walks 
              S(i,:)=best+0.001*randn(1,d);
          end

     % Evaluate new solutions
          S(i,:)=simplebounds(S(i,:),Lb,Ub);
           Fnew=shannonEntropy(S(i,:),h);
     % Update if the solution improves, or not too loud
           if (Fnew>=Fitness(i)) & (rand<A) ,
                Sol(i,:)=S(i,:);
                Fitness(i)=Fnew;
           end

          % Update the current best solution
          if Fnew>=fmin,
                best=S(i,:);
                fmin=Fnew;
          end
        end
        N_iter=N_iter+n;
end
best=simplebounds(best,Lb,Ub);
fitness(pp)=fmin;
end;
mea=mean(fitness);
standereddeviation=std(fitness);
maxfitness=max(fitness);
% Output/display
disp(['Number of evaluations: ',num2str(N_iter)]);
disp(['Best =',num2str(best),' fmin=',num2str(fmin)]);

% Application of simple limits/bounds
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound vector
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bound vector 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;
  
% function entropy=shannonEntropy(x,h)
% x=round(x);
% entropy=0;
% x=[1 (x+1) 256];
% % x=256.*ones(1,6);
% x(x>256)=256;
% s=size(x);
% for kk=1:size(x,2)
%     if x(kk)<0
%         x(kk)=-x(kk);
%     end;
% end;
% x=sort(x,'ascend');
% for i=1:s(2)-1
%     x(x>256)=256;
%     temp = h(x(i):x(i+1));
%     tsum= sum(temp);
%     ent=0;
%     if tsum~=0
%     for j=x(i):x(i+1)
%         if h(j)~=0
%             a=h(j)/tsum;
%             ent = ent + a*log(a);
%         end
%     end
%     end
%     entropy = entropy + ent;    
% end
% 
% entropy=-entropy;
