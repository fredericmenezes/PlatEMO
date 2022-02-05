%--------------------------------------------------------------------------%
% This function can be referred in the below book:                     %
% Levy exponent and coefficient                                            %
% For details, see Chapter 3 of the following book:                        %
% Xin-She Yang, Nature-Inspired Optimization Algorithms, Elsevier, (2014). %
% Chaper 3                                                                 %
%--------------------------------------------------------------------------%
function s=Levy(d) % d is number of dimension
beta=3/2;
% Eq. (3.27)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
% Eq. (3.26)
u=randn(1,d)*sigma;
v=randn(1,d);
% Eq. (3.25)
step=u./abs(v).^(1/beta);
s=0.01*step;