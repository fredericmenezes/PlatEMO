function POS = opposJumping(POS, var_min, var_max, it, MaxIt)
% Update new position by opposition-based jumping using Eq. 26
%
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if rand <(it/MaxIt)
        [POS] = OPPOS(POS,var_max,var_min);
    end
    %% Check boundaries
    if rand>=0.5
        POS=checksimplebounds(POS,var_min',var_max');
    else
        POS = checkBoundaries(POS,var_max,var_min);
    end
end

%% This function checks the boundary of jellyfish search space
function [POS] = checkBoundaries(POS,var_max,var_min)
    % Useful matrices
    Np = size(POS,1);
    MAXLIM   = repmat(var_max(:)',Np,1);
    MINLIM   = repmat(var_min(:)',Np,1);
    POS(POS>MAXLIM) = MAXLIM(POS>MAXLIM);
    POS(POS<MINLIM) = MINLIM(POS<MINLIM);
end
function POS=checksimplebounds(POS,Lb,Ub)
    for i=1:size(POS,1)
        ns_tmp=POS(i,:);
        I=ns_tmp<Lb;
        while sum(I)~=0
            ns_tmp(I)=Ub(I)+(ns_tmp(I)-Lb(I));
            I=ns_tmp<Lb;
        end
        J=ns_tmp>Ub;
        while sum(J)~=0
            ns_tmp(J)=Lb(J)+(ns_tmp(J)-Ub(J));
            J=ns_tmp>Ub;
        end
        POS(i,:)=ns_tmp;
    end
end