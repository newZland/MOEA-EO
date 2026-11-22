function Offspring = OP_SBX(Parent,Args)
% SBX operator.

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

    [prob,dist] = deal(Args{:});
    [N,D] = size(Parent);
    Offspring = Parent;
    if N ~= 0 && D ~= 0
        Parent_ref = [Parent(end,:);Parent(1:end-1,:)];
        beta = zeros(N,D);
        mu = rand(N,D);
        beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(dist+1));
        beta(mu>0.5) = (2-2*mu(mu>0.5)).^(-1/(dist+1));
        beta = beta.*(-1).^randi([0,1],N,D);
        beta(rand(N,D)<0.5) = 1;
        beta(repmat(rand(N,1)>prob,1,D)) = 1;
        k = rand(N,D) <= 0.5;
        Offspring(k==1) = (Parent(k==1)+Parent_ref(k==1))/2+beta(k==1).*(Parent(k==1)-Parent_ref(k==1))/2;
        Offspring(k~=1) = (Parent(k~=1)+Parent_ref(k~=1))/2-beta(k~=1).*(Parent(k~=1)-Parent_ref(k~=1))/2;
    end
end