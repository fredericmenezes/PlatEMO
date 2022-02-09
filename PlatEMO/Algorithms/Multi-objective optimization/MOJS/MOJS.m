%classdef MOJS < ALGORITHM
% <multi> <real>
% Multi-Objective Jellyfish Search
% ngrid ---  20 ---  Numero de grades em cada dimensão
% Nr    --- 100 ---  Tamanho do Arquivo
% MaxIt ---  99 ---  Maximo de iteracoes. Se alterar maxFE faça MaxIt = (maxFE-N)/N


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
            % Parameters
            Np      = Problem.N;
            [ngrid,Nr, MaxIt] = Algorithm.ParameterSet(20,100,99);
            nVar    = Problem.D;
            var_min = Problem.lower(:);
            var_max = Problem.upper(:);
            it=1;
            % Initialization by Eq. 25
            Population = Problem.InitializationChaos();

            POS          = Population.decs;
            POS_fit      = Population.objs;
            ELI_POS      = POS;
            ELI_POS_fit  = POS_fit;
            DOMINATED    = checkDomination(POS_fit);

            Archive      = Population(~DOMINATED);
            ARCH.pos     = Archive.decs;
            ARCH.pos_fit = Archive.objs;
            ARCH         = updateGrid(ARCH,ngrid);
            
            %display(['Iteration #0 - Archive size: ' num2str(size(ARCH.pos,1))]);
            %% Main MOJS loop
            while Algorithm.NotTerminated(Archive)

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
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Update new position by opposition-based jumping using Eq. 26
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if rand <(it/MaxIt)
                    [POS] = OPPOS(POS,var_max,var_min);
                end
                %% Check boundaries
                if rand>=0.5
                    POS=checksimplebounds(POS,var_min',var_max');
                else
                    POS = checkBoundaries(POS,var_max,var_min);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Evaluate the population and update the archive
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                %display(['Iteration #' num2str(it) ' - Evaluations #' num2str(Problem.FE) ' - Archive size: ' num2str(size(ARCH.pos,1))]);
                
                Archive = SOLUTION(ARCH.pos);
                it=it+1;

            end

           
        end
    end
end