%% This function updates the archive given a new population 
function ARCH = updateArchive(ARCH,POS,POS_fit,ngrid)
    % Domination between jellyfish
    DOMINATED  = checkDomination(POS_fit);
    ARCH.pos    = [ARCH.pos; POS(~DOMINATED,:)];
    ARCH.pos_fit= [ARCH.pos_fit; POS_fit(~DOMINATED,:)];
    % Domination between nondominated jellyfish and the last archive 
    DOMINATED  = checkDomination(ARCH.pos_fit);
    ARCH.pos_fit= ARCH.pos_fit(~DOMINATED,:);
    ARCH.pos    = ARCH.pos(~DOMINATED,:);
    % Updating the grid
    ARCH        = updateGrid(ARCH,ngrid);
end