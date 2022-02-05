classdef MOJS < ALGORITHM
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

	methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [ngrid,Nr] = Algorithm.ParameterSet(20,100);
            var_min    = Problem.lower(:);
            var_max    = Problem.upper(:);
            MaxIt      = ceil((Problem.maxFE-Problem.N)/Problem.N);
            if MaxIt<1; MaxIt=1; end
            
            %% Generate random population
            %  Initialization by Eq. 25
            Population = Problem.InitializationChaos();
            
            POS          = Population.decs;
            POS_fit      = Population.objs;
            ELI_POS      = POS;
            ELI_POS_fit  = POS_fit;
            DOMINATED    = checkDomination(POS_fit);
            ARCH.pos     = POS(~DOMINATED,:);
            ARCH.pos_fit = POS_fit(~DOMINATED,:);
            ARCH         = updateGrid(ARCH,ngrid);
            
            it=1;
            %% Optimization
            while Algorithm.NotTerminated(Archive)
                POS = OperatorJS(POS,ELI_POS,ELI_POS_fit, ARCH,it,MaxIt);
                POS = opposJumping(POS,var_min,var_max,it,MaxIt);
                [ELI_POS, ELI_POS_fit] = evaluatePOS(POS,ELI_POS,ELI_POS_fit);
                
                
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
                
                it=it+1;
                
            end
        end
    end
end