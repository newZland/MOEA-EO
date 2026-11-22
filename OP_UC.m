% REAL类型，均匀交叉，子代随机从父代集成某一特征
% Parent父类，prob交叉概率
function Offspring = OP_UC(Parent,Args)
% The environmental selection of MOEA-EO.

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

    [prob] = deal(Args{:});
    [N,D] = size(Parent);
    Offspring = Parent;
    if N ~= 0 && D ~= 0
        Parent_ref = [Parent(end,:);Parent(1:end-1,:)];
        k = rand(N,D) < prob;
        Offspring(k) = Parent_ref(k);
    end
end