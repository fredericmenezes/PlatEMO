function selected = selectLeader(ARCH)
% This function calucates the leader performance by a roulette wheel 
% selection based on the quality of each hypercube

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% Roulette wheel
prob    = cumsum(ARCH.quality(:,2));     % Cumulated probs
sel_hyp = ARCH.quality(find(rand(1,1)*max(prob)<=prob,1,'first'),1); % Selected hypercube
% Select the index leader as a random selection inside that hypercube
idx      = 1:1:length(ARCH.grid_idx);
selected = idx(ARCH.grid_idx==sel_hyp);
selected = selected(randi(length(selected)));
end