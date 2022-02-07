function [ELI_POS, ELI_POS_fit] = evaluatePOS(POS, ELI_POS, ELI_POS_fit, Problem)
% Evaluate the population
%
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    Np = size(POS,1);
    
    POS_fit = Problem.CalObj(POS);
    
    pos_best = dominates(POS_fit, ELI_POS_fit);
    best_pos = ~dominates(ELI_POS_fit, POS_fit);
    
    best_pos(rand(Np,1)>=0.5) = 0;
    
    if (sum(pos_best)>1)
        ELI_POS_fit(pos_best,:) = POS_fit(pos_best,:);
        ELI_POS(pos_best,:) = POS(pos_best,:);
    end
    if (sum(best_pos)>1)
        ELI_POS_fit(best_pos,:) = POS_fit(best_pos,:);
        ELI_POS(best_pos,:) = POS(best_pos,:);
    end
    
end