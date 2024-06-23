% Interior search algorithm (ISA)
function [Thresholds,meann,stdd,maxfitness]=ISA(h,Level,N_iter)
nr=30; % Number of runs
n = 10*Level; % Population size
%N_iter = N_iter; % Number of iterations

% Lower and upper bounds
Lb = 0;
Ub = 255;

d = (Level); % Problem dimesion

% Initial Matrices
ns=zeros(n,d); x_new=ns; Cmax=0; Cn = zeros(1,2*n);
fvalue=zeros(1,n); fvalue_new=zeros(1,n); Xibest=zeros(d,1);
f_history=zeros(1,N_iter); RunHistory=zeros(N_iter,nr);
for r = 1:nr % number of runs
    % Random initial solutions and evaluation of them
    for i=1:n
        ns(i,:) = Lb+rand(1,d).*(Ub-Lb);
        fvalue(i) = shannonEntropy(ns(i,:),h');
        C(i,:) = constraints(ns(i,:));
    end

    % Iterations begin
    for j=1:N_iter,
        alpha = (j/N_iter);
        LX=min(ns);
        UX=max(ns);
        [~,l] = max(fvalue);
        center = ns(l,:);
        for i = 1:n
            if i == l
                x_new(i,:) = ns(i,:)+randn(1,d).*(Ub-Lb).*0.01;
            else
                if rand < alpha
                    beta = rand;
                    mirror = beta*ns(i,:)+(1-beta)*center;
                    x_new(i,:) = 2*mirror - ns(i,:);
                else
                    x_new(i,:) = LX+(UX-LX).*rand(1,d);
                end
            end

            % % % % Apply the evolutionary boundry handling scheme % % % % 
            ns_tmp=x_new(i,:);
            I=ns_tmp<Lb;
            A=rand;
            ns_tmp(I)=A*Lb+(1-A)*center(I);
            J=ns_tmp>Ub;
            B=rand;
            ns_tmp(J)=B*Ub+(1-B)*center(J);
            x_new(i,:)=ns_tmp;
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            fvalue_new(i)=shannonEntropy(x_new(i,:),h');
            C(n+i,:) = constraints(x_new(i,:));
        end

        % % Apply the parameter-less penalty constraint handling scheme % %
        C = (C >= 0).*C;
        Fmax = max(max(fvalue,fvalue_new));
        Cmax = max(max(C),Cmax);
        
        for i=1:n
            CC=C(i,:)./Cmax; CC(isnan(CC))=0; Cn(i) = sum(CC);
            if 0 < Cn(i)
                fvalue(i) = Fmax + Cn(i);
            end
            CC_new=C(n+i,:)./Cmax; CC_new(isnan(CC_new))=0; Cn(n+i) = sum(CC_new);
            if 0 < Cn(n+i)
                fvalue_new(i) = Fmax + Cn(n+i);
            end
            if fvalue_new(i)>fvalue(i)
                ns(i,:)=x_new(i,:);
                fvalue(i)=fvalue_new(i);
            end
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        % storing the statistical analysis of fvalue in each iteration
        f_history(j)= max(fvalue);

    end % end of iterations
    % Post-processing and storing the runs history
    [~, l] = max(fvalue);
    Xibest(:,r) = ns(l,:);
    RunHistory(:,r) = f_history;
     [Bestf Bestrun] = max(RunHistory(end,:));
     RunHistoryy(nr)=Bestf;
end
%% Post Processing
[Bestf Bestrun] = max(RunHistory(end,:))
Thresholds = round(Xibest(:,Bestrun))
meann =  mean(RunHistoryy)
Medianf =  median(RunHistoryy)
maxfitness =  max(RunHistoryy)
stdd = std(RunHistoryy)



% Nonlinear constraints
function c=constraints(x)
format long
% Spring
c(1)=1-x(2)^3*x(3)/(71785*x(1)^4);
c(2)=(4*x(2)^2-x(1)*x(2))/(12566*(x(2)*x(1)^3-x(1)^4))+1/(5108*x(1)^2)-1;
c(3)=1-140.45*x(1)/(x(2)^2*x(3));
c(4)=(x(1)+x(2))/1.5-1;
if isempty(c);c = -1;end % For an unconstraint function