function [T,meann,stdd,minimum]=pso(h,D,max_runs)
range=[1 250];
n_var=D; %Dimensions
numInd=10*D;  %%%population
pesoStoc=0.1;
numIter =max_runs;
range_min=range(1); % Range for initial swarm's elements
range_max=range(2);
numVar=n_var; % Number of variables
%fOb=func; % Objective function
for pp=1:max_runs
iter=1; % Number of iteration
ind = range_min + (range_max-range_min).*rand(numInd,n_var); % Initial swarm
k=pesoStoc; % weight of stocastic element

v=zeros(numInd,numVar); % Vector of swarm's velocity

radius=1000; % Initial radius for stop criterion
tolerance=1;
while iter<=numIter% && radius>tolerance
    for l=1:numInd
        %valF(l,1)=fOb(ind(l,:)); % Fitness function for the swarm
         ind(l,:)=sort(ind(l,:));
        valF(l,1)=fuzzyEntropy(ind(l,:),h');
    end
    [valF_ord,index]=sort(valF,'descend'); % Sort the objective function's values for the swarm and identify the leader
    leader=ind(index(1),:);
    fmin(iter)=max(valF);
    for l=1:size(ind,1) % Calculates the new velocity and positions for all swarm's elements
        fi=rand();
        v(l,:)=(1-(sqrt(k*fi))/2)*v(l,:)+k*fi*(leader-ind(l,:)); % Velocity
        ind(l,:)=ind(l,:)+(1-(sqrt(k*fi))/2)*v(l,:)+(1-k*fi)*(leader-ind(l,:)); % Position
         
    end
    radius=norm(leader-ind(index(end),:)); % Calculates the new radius
    iter=iter+1; % Increases the number of iteration
   ind(ind<0)=50.*rand(1,1);
        ind(ind>255)=255.*rand(1,1);
end
p_min=valF; % Output variables
f_min(pp)=max(fmin);

end
minimum=max(f_min)
meann=mean(f_min)
stdd=std(f_min)
T=leader;