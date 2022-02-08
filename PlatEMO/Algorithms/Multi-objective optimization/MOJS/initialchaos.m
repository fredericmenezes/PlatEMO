%------------------------------------------------------------------------------------------------------------------------%
% The equations can be referred to the below paper:                                                                      %
% Fister I, Perc M, Kamal SM, Fister I. A review of chaos-based firefly algorithms: Perspectives and research challenges % 
% Applied Mathematics and Computation                                                                                    %
% 2015;252:155-65.DOI: https://doi.org/10.1016/j.amc.2014.12.006                                                         %
%------------------------------------------------------------------------------------------------------------------------%
function pop=initialchaos(index,num_pop,nd,Ub,Lb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the population by chaostic map
% This program is written at PiM lab
% pop=initialchaos(index,num_pop,nd,Ub,Lb)
% index: choose the map case
%       1:Chebyshev map
%       2:Circle map
%       3:Gauss/mouse map
%       4:Intermittency map
%       5:Iterative map
%       6:Liebovitch map
%       7:Logistic map
%       8:Piecewise map
%       9:Sine map
%       10:Singer map
%       11:Sinusoidal map
%       12:Tent map
%       13: Kent map
% num_iter: Number of population;
% nd: Number of dimention; e.g: nd=4;
% Ub: Matrix of Upper bound,e.g:[1 1 1 1];
% Lb: Matrix of lower bound,e.g:[-1 0 -2 3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initinationthearms(n,d,Ub,Lb)
%chaos(chaosIndex-1,iteration,max_it,chValue);
if size(Lb,2)==1
    Lb=Lb*ones(1,nd);
    Ub=Ub*ones(1,nd);
end
x(1,:)=rand(1,nd);%0.7;
switch index
%% Chebyshev map
    case 1
        for i=1:(num_pop-1)
            x(i+1,:)=cos(i*acos(x(i,:)));
        end
%% Circle map
    case 2
        a=0.5;
        b=0.2;
        for i=1:(num_pop-1)
            x(i+1,:)=mod(x(i,:)+b-(a/(2*pi))*sin(2*pi*x(i,:)),1);
        end
%% Gauss/mouse map
    case 3
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)==0
                    x(i+1,k)=0;
                else
                    x(i+1,k)=mod(1/x(i,k),1);
                end
            end
        end
%% Intermittency map
    case 4
        P=0.4;
        m=2;
        eps=1e-10;
        c=(1-eps-P)/(P^m);
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)<=P
                    x(i+1,k)=eps+x(i,k)+c*(x(i,k))^m;
                else
                    x(i+1,k)=(x(i,k)-P)/(1-P);
                end
            end
        end
%% Iterative map
    case 5   
        a=0.7;
        for k=1:nd
            for i=1:(num_pop-1)
                x(i+1,k)=sin((a*pi)/x(i,k));
            end
        end
%% Liebovitch map 
    case 6
        
        P1=0.3;
        P2=0.7;
        anpha=P2/P1*(1-(P2-P1));
        beta=(1/(P2-1))*((P2-1)-P1*(P2-P1));
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)<=P1
                    x(i+1,k)=anpha*x(i,k);
                elseif x(i,k)<=P2
                    x(i+1,k)=(P2-x(i,k))/(P2-P1);
                else
                    x(i+1,k)=1-beta*(1-x(i,k));
                end
            end
        end
%% Logistic map
    case 7 
        a=4;
        for i=1:(num_pop-1)
            x(i+1,:)=a*x(i,:).*(1-x(i,:));
        end
%% Piecewise map
    case 8   
        P=0.4;
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)>=0 && x(i,k)<P
                    x(i+1,k)=x(i,k)/P;
                end
                if x(i,k)>=P && x(i,k)<0.5
                    x(i+1,k)=(x(i,k)-P)/(0.5-P);
                end
                if x(i,k)>=0.5 && x(i,k)<1-P
                    x(i+1,k)=(1-P-x(i,k))/(0.5-P);
                end
                if x(i,k)>=1-P && x(i,k)<1
                    x(i+1,k)=(1-x(i,k))/P;
                end
            end
        end
%% Sine map
    case 9      
        for i=1:(num_pop-1)
            x(i+1,:) = sin(pi*x(i,:));
        end
%% Singer map
    case 10       
        u=1.07;
        for k=1:nd
            for i=1:(num_pop-1)
                x(i+1,k) = u*(7.86*x(i,k)-23.31*(x(i,k)^2)+28.75*(x(i,k)^3)-13.302875*(x(i,k)^4));
            end
        end
%% Sinusoidal map
    case 11        
        for k=1:nd
            for i=1:(num_pop-1)
                while mod(x(i,k),1)==0
                    x(i,k)=rand(1,1);
                end
                x(i+1,k) = 2.3*x(i,k)^2*sin(pi*x(i,k));
            end
        end
%% Tent map
    case 12        
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)<0.7
                    x(i+1,k)=x(i,k)/0.7;
                end
                if x(i,k)>=0.7
                    x(i+1,k)=(10/3)*(1-x(i,k));
                end
            end
        end
%% Kent map
    case 13        
        m=0.6;
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)<=m
                    x(i+1,k)=x(i,k)/m;
                else
                    x(i+1,k)=(1-x(i,k))/(1-m);
                end
            end
        end
end
for k=1:nd
    for i=1:num_pop
        pop(i,k)=Lb(k)+x(i,k)*(Ub(k)-Lb(k));
    end
end
end