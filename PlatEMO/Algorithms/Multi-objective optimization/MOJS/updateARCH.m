function Archive = updateARCH(ARCH, POS, POS_fit, ELI_POS, ELI_POS_fit, ngrid, Nr)
% Update the archive  
% 
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

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
    Archive = SOLUTION(ARCH.pos,ARCH);
end

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
%% This function deletes an excess of jellyfish inside the archive using crowding distances
function ARCH = deleteFromArchive(ARCH,n_extra,ngrid)
    % Compute the crowding distances
    crowding = zeros(size(ARCH.pos,1),1);
    for m = 1:1:size(ARCH.pos_fit,2)
        [m_fit,idx] = sort(ARCH.pos_fit(:,m),'ascend');
        m_up     = [m_fit(2:end); Inf];
        m_down   = [Inf; m_fit(1:end-1)];
        distance = (m_up-m_down)./(max(m_fit)-min(m_fit));
        [~,idx]  = sort(idx,'ascend');
        crowding = crowding + distance(idx);
    end
    crowding(isnan(crowding)) = Inf;
    % This function deletes the extra jellyfish with the smallest crowding distances
    [~,del_idx] = sort(crowding,'ascend');
    del_idx = del_idx(1:n_extra);
    ARCH.pos(del_idx,:) = [];
    ARCH.pos_fit(del_idx,:) = [];
    ARCH = updateGrid(ARCH,ngrid);
end
