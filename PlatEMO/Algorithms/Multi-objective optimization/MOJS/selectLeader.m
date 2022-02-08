%% This function calucates the leader performance by a roulette wheel selection
% based on the quality of each hypercube
function selected = selectLeader(ARCH)
    % Roulette wheel
    prob    = cumsum(ARCH.quality(:,2));     % Cumulated probs
    sel_hyp = ARCH.quality(find(rand(1,1)*max(prob)<=prob,1,'first'),1); % Selected hypercube
    % Select the index leader as a random selection inside that hypercube
    idx      = 1:1:length(ARCH.grid_idx);
    selected = idx(ARCH.grid_idx==sel_hyp);
    selected = selected(randi(length(selected)));
end