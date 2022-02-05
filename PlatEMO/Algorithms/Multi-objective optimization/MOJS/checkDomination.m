function domi_vector = checkDomination(fitness)
% This function checks the domination inside the population.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    Np = size(fitness,1);
    if Np>2
        domi_vector = zeros(Np,1);
        all_perm = nchoosek(1:Np,2);    % Possible permutations
        all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];

        d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
        dominated_particles = unique(all_perm(d==1,2));
        domi_vector(dominated_particles) = 1;
    else
        domi_vector=ones(Np,1);
    end
end

