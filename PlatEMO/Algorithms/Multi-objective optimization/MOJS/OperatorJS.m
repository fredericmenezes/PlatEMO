function POS = OperatorJS(POS, ELI_POS, ELI_POS_fit, ARCH, it, MaxIt)
%OperatorPSO - The operator of jellyfish search optimization.
%
% POS = OperatorJS(POS, ELI_POS, ELI_POS_fit, ARCH, it, MaxIt) uses the 
% operator of jellyfish search optimization to generate population-based 
% on population POS, elite population ELI_POS, objective values of elite 
% population ELI_POS_fit, and an archive of the best jellyfishes ARCH. POS 
% and ELI_POS are matrices, and ARCH should be a struct that contains:
%
%       ARCH.pos
%       ARCH.pos_fit
%       ARCH.hypercube_limits
%       ARCH.grid_subidx
%       ARCH.quality
%
%   Example:
%       POS = OperatorPSO(POS,ELI_POS, ELI_POS_fit,ARCH,it,MaxIt)

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

[Np, nVar] = size(POS);
% Select leader by Eq. 16
h = selectLeader(ARCH);
% Calculate time control by Eq. 15
Ct=abs((1-it*((1)/MaxIt))*(2*rand-1));
if Ct>=0.5
    Meanvl=mean(ELI_POS);
    for i=1:Np
        % The new position is determined by Eq.19 and Eq.20
        POS(i,:) = ELI_POS(i,:) + Levy(nVar).*(ARCH.pos(h,:) - 3*rand([1 nVar]).*Meanvl);
    end
else
    for i=1:Np
        if rand<(1-Ct)
            % Jellyfish follow type B
            % Determine the direction by Eq. 24
            j=i;
            while j==i
                j=randperm(Np,1);
            end
            Step = ELI_POS(i,:) - ELI_POS(j,:);
            if dominates(ELI_POS_fit(j,:),ELI_POS_fit(i,:))
                Step = -Step;
            end
            % The new position is determined by Eq. 22
            POS(i,:) =ARCH.pos(h,:) + rand([1 nVar]).*Step;
        else
            % Jellyfish follow type A
            % The new position is determined by Eq. 21
            POS(i,:)=ARCH.pos(h,:)+Levy(nVar).*(ELI_POS(i,:)-ARCH.pos(h,:));
        end
    end
end

end