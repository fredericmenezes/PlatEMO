%----------------------------------------------------------------------------------------%
% The equations can be referred to the below paper:                                      %
% Equation 01.                                                                           %
% Paper: Mahdavi S, Rahnamayan S, Deb K. Opposition based learning: A literature review. % 
% Swarm and Evolutionary Computation 2018;39:1-23.                                       %
% DOI:  https://doi.org/10.1016/j.swevo.2017.09.010                                      %
%----------------------------------------------------------------------------------------% 

function [POS] = OPPOS(POS,var_max,var_min)
    Np = size(POS,1);
    MAXLIM   = repmat(var_max(:)',Np,1);
    MINLIM   = repmat(var_min(:)',Np,1);
    POS = (MINLIM+MAXLIM)-POS;
end