%classdef MOJS < ALGORITHM
% <multi> <real>
% Multi-Objective Jellyfish Search
% ngrid ---  20 ---  Number of grids in each dimension
% Nr    --- 100 ---  Archive size

%------------------------------- Reference --------------------------------
% J. S. Chou, D. N. Truong, Multiobjective optimization inspired by 
% behavior of jellyfish for solving structural design problems, Chaos,
% Solitons & Fractals, 2020, 135, 109738.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
classdef MOJS < ALGORITHM
    methods
        function main(Algorithm,Problem)
            
            %g = @(x) 1 + 9*mean(x(:,2:end),2);
            %h = @(x) 1 - (x(:,1)./g(x)).^0.5;
            %fun = @(x) [x(:,1), g(x).*h(x)];
            
            
            % Parameters
            Np      = Problem.N;
            %Nr      = params.Nr;
            %MaxIt   = 400;
            [ngrid,Nr, MaxIt] = Algorithm.ParameterSet(20,100,400);
            %ngrid   = params.ngrid;
            %fun     = MultiObj.fun;
            nVar    = Problem.D;
            var_min = Problem.lower(:);
            var_max = Problem.upper(:);
            it=1;
            % Initialization by Eq. 25
            POS=initialchaos(7,Np,nVar,var_max',var_min');
            POS_fit      = Problem.CalObj(POS);
            ELI_POS      = POS;
            ELI_POS_fit  = POS_fit;
            DOMINATED    = checkDomination(POS_fit);
            ARCH.pos     = POS(~DOMINATED,:);
            ARCH.pos_fit = POS_fit(~DOMINATED,:);
            ARCH         = updateGrid(ARCH,ngrid);
            Archive      = SOLUTION(ARCH.pos);
            
            display(['Iteration #0 - Archive size: ' num2str(size(ARCH.pos,1))]);
            %% Main MOJS loop
            stopCondition = false;
            while (~stopCondition || Algorithm.NotTerminated(Archive))
                % Select leader by Eq. 16
                h = selectLeader(ARCH);
                % Calculate time control by Eq. 15
                Ct=abs((1-it*((1)/MaxIt))*(2*rand-1));
                if Ct>=0.5
                    Meanvl=mean(ELI_POS);
                    for i=1:Np
                        % The new position is determined by Eq.19 and Eq.20
                        POS(i,:) = ELI_POS(i,:) + Levy(nVar).*(ARCH.pos(h,:) - 3*rand([1 nVar]).*Meanvl);
                    end
                else
                    for i=1:Np
                        if rand<(1-Ct)
                            % Jellyfish follow type B
                            % Determine the direction by Eq. 24
                            j=i;
                            while j==i
                                j=randperm(Np,1);
                            end
                            Step = ELI_POS(i,:) - ELI_POS(j,:);
                            if dominates(ELI_POS_fit(j,:),ELI_POS_fit(i,:))
                                Step = -Step;
                            end
                            % The new position is determined by Eq. 22
                            POS(i,:) =ARCH.pos(h,:) + rand([1 nVar]).*Step;
                        else
                            % Jellyfish follow type A
                            % The new position is determined by Eq. 21
                            POS(i,:)=ARCH.pos(h,:)+Levy(nVar).*(ELI_POS(i,:)-ARCH.pos(h,:));
                        end
                    end
                end
                %% Update new position by opposition-based jumping using Eq. 26
                if rand <(it/MaxIt)
                    [POS] = OPPOS(POS,var_max,var_min);
                end
                %% Check boundaries
                if rand>=0.5
                    POS=checksimplebounds(POS,var_min',var_max');
                else
                    POS = checkBoundaries(POS,var_max,var_min);
                end
                %% Evaluate the population
                POS_fit = Problem.CalObj(POS);
                pos_best = dominates(POS_fit, ELI_POS_fit);
                best_pos = ~dominates(ELI_POS_fit, POS_fit);
                best_pos(rand(Np,1)>=0.5) = 0;
                if(sum(pos_best)>1)
                    ELI_POS_fit(pos_best,:) = POS_fit(pos_best,:);
                    ELI_POS(pos_best,:) = POS(pos_best,:);
                end
                if(sum(best_pos)>1)
                    ELI_POS_fit(best_pos,:) = POS_fit(best_pos,:);
                    ELI_POS(best_pos,:) = POS(best_pos,:);
                end
                %% Update the archive
                if size(ARCH.pos,1)==1
                    ARCH.pos= POS;
                    ARCH.pos_fit= POS_fit;
                    ARCH = updateArchive(ARCH,ELI_POS,ELI_POS_fit,ngrid);
                else
                    ARCH = updateArchive(ARCH,ELI_POS,ELI_POS_fit,ngrid);
                    if size(ARCH.pos,1)==1
                        ARCH.pos= ELI_POS;
                        ARCH.pos_fit= ELI_POS_fit;
                    end
                end
                if(size(ARCH.pos,1)>Nr)
                    % Delete the worst members from archive by Eq. 18
                    ARCH = deleteFromArchive(ARCH,size(ARCH.pos,1)-Nr,ngrid);
                end
                display(['Iteration #' num2str(it) ' - Evaluations #' num2str(Problem.FE) ' - Archive size: ' num2str(size(ARCH.pos,1))]);
                
                Archive = SOLUTION(ARCH.pos);
                it=it+1;
                if(it>MaxIt), stopCondition = true; end
            end
            %% Plotting paretofront
            
            %figure
            
            %ObjectiveFunction=fun;
            %x=MultiObj.var_min(1):0.01:MultiObj.var_max(1);
            
            %for i=1:size(x,2)
            %    TPF(i,:)=ObjectiveFunction(x(i));
            %end
            
            %line(TPF(:,1),TPF(:,2));
            %title('Schaffer')
            
            %xlabel('f1')
            %ylabel('f2')
            
            %box on
            
            %fig=gcf;
            
            %set(findall(fig,'-property','FontName'),'FontName','Garamond')
            %set(findall(fig,'-property','FontAngle'),'FontAngle','italic')
            
            %hold on
            
            if(size(ARCH.pos_fit,2)==2)
                plot(ARCH.pos_fit(:,1),ARCH.pos_fit(:,2),'or'); hold on;
                grid on; xlabel('f1'); ylabel('f2'); %legend('PF verdadeiro','MOJS');
            end
            
            if(size(ARCH.pos_fit,2)==3)
                plot3(ARCH.pos_fit(:,1),ARCH.pos_fit(:,2),ARCH.pos_fit(:,3),'or'); hold on;
                grid on; xlabel('f1'); ylabel('f2'); zlabel('f3');
            end
            %set(gca, 'XLim', [0 4.5])
            %set(gca, 'YLim', [0 4.5])
            disp('fred');
            
            
            
        end
    end
end