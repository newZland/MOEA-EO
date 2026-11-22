function Offspring = OP_DE(Parent,Args)
% DE operator.

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

    [prob,ratio] = deal(Args{:});
    [N,D] = size(Parent);
    Offspring = Parent;
    if N ~= 0 && D ~= 0
        de1 = Parent(randperm(N),:);
        de2 = Parent(randperm(N),:);
        sle = rand(N,D) < prob;
        Offspring(sle) = Offspring(sle) + ratio*(de1(sle)-de2(sle));
    end
end