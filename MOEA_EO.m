classdef MOEA_EO < ALGORITHM
% <2026> <multi/many> <real> <large/none> <sparse/none>
% Multi-objective evolutionary algorithm with nonzero detection

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

    methods
        function main(Algorithm,Problem)
            % Generate the initial population and reference points, and initialize the ensemble operator handler
            Operators = {{'OP_SBX',{1,0.5}},{'OP_UC',{0.5}},{'OP_SDM',{0.01}},{'OP_PM',{1,5}},{'OP_DE',{1,1.1}}};
            EO = HC(Problem,Operators);
            Population = InitPopulation(Problem,0.75,1);
            [Z,~] = UniformPoint(Problem.N,Problem.M);
            Zmin = min(Population(all(Population.cons<=0,2)).objs,[],1);
            % Optimization
            while Algorithm.NotTerminated(Population)
               Offspring  = EO.Generate(Population,Problem);
               Zmin = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
               Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
            end
        end
    end
end

