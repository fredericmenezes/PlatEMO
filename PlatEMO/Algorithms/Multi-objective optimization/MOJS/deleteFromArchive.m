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
