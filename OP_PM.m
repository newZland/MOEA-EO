function Offspring = OP_PM(Parent,Args)
% PM operator.

%------------------------------- Reference --------------------------------
% D. Liu, X. Li, A reproduction operator ensemble evolutionary algorithm 
% for sparse large-scale multi-objective optimization, Applied Soft 
% Computing, 2026, 1–25, doi: https://doi.org/10.1016/j.asoc.2026.116884.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Dun Liu (email: liudun@zjnu.edu.cn)

    [lower,upper,prob,dist] = deal(Args{:});
    [N,D] = size(Parent);
    Offspring = Parent;

    if N ~= 0 && D ~= 0
        lower = repmat(lower,N,1);
        upper = repmat(upper,N,1);
        sle = rand(N,D) < prob/D;
        mu = rand(N,D);
        ch  = sle & mu<=0.5;
        Offspring = min(max(Parent,lower),upper);
        Offspring(ch) = Offspring(ch)+(upper(ch)-lower(ch)).*((2.*mu(ch)+(1-2.*mu(ch)).*(1-(Offspring(ch)-lower(ch))./(upper(ch)-lower(ch))).^(dist+1)).^(1/(dist+1))-1);
        ch = sle & mu>0.5; 
        Offspring(ch) = Offspring(ch)+(upper(ch)-lower(ch)).*(1-(2.*(1-mu(ch))+2.*(mu(ch)-0.5).*(1-(upper(ch)-Offspring(ch))./(upper(ch)-lower(ch))).^(dist+1)).^(1/(dist+1)));
    end
end