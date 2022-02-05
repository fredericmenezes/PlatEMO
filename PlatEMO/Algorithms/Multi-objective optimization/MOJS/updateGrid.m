function ARCH = updateGrid(ARCH,ngrid)
% Function that updates the hypercube grid, the hypercube where belongs

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Computing the  hypercube limitation
    ndim = size(ARCH.pos_fit,2);
    ARCH.hypercube_limits = zeros(ngrid+1,ndim);
    for dim = 1:1:ndim
        ARCH.hypercube_limits(:,dim) = linspace(min(ARCH.pos_fit(:,dim)),max(ARCH.pos_fit(:,dim)),ngrid+1)';
    end
    % Computing where belongs each jellyfish
    npar = size(ARCH.pos_fit,1);
    ARCH.grid_idx = zeros(npar,1);
    ARCH.grid_subidx = zeros(npar,ndim);
    for n = 1:1:npar
        idnames = [];
        for d = 1:1:ndim
            ARCH.grid_subidx(n,d) = find(ARCH.pos_fit(n,d)<=ARCH.hypercube_limits(:,d)',1,'first')-1;
            if(ARCH.grid_subidx(n,d)==0), ARCH.grid_subidx(n,d) = 1; end
            idnames = [idnames ',' num2str(ARCH.grid_subidx(n,d))];
        end
        ARCH.grid_idx(n) = eval(['sub2ind(ngrid.*ones(1,ndim)' idnames ');']);
    end
    % Quality based on the number of jellyfish in each hypercube
    ARCH.quality = zeros(ngrid,2);
    ids = unique(ARCH.grid_idx);
    for i = 1:length(ids)
        ARCH.quality(i,1) = ids(i);                       
        ARCH.quality(i,2) = 10/sum(ARCH.grid_idx==ids(i)); 
    end
end

