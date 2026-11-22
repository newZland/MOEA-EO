classdef HC < handle
%HC - The HC structure of MOEA-EO.
%   This is a superclass for the horizontal concatenation structure of
%   MOEA-EO.
%
% HC properties:
%   op           <scalar>   operator pool
%   cnt          <scalar>   number of parents assigned to each operator
%
% HC methods:
%   HC           <public>   the constructor setting all the properties specified by user
%   Generate     <public>   generate offspring from parents

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

% This function is rewritten by Dun Liu (email: liudun@zjnu.edu.cn)

    properties    
        op;         % Operator pool
        cnt;      % Number of parents assigned to each operator
    end
    methods
        function obj = HC(prob,OP)
        %HC - Initialize a handle to the HC structure given the operator pool OP
        %
        %   Example:
        %       EO = HC(Problem,{{'OP_SBX',{1,0.5}},{'OP_UC',{0.5}},{'OP_SDM',{0.01}},{'OP_PM',{1,5}},{'OP_DE',{1,0.9}}});
            len = size(OP,2);
            if len ~= 0
                obj.op = OP;
                obj.cnt = floor(prob.N/len)*ones(1,len);
                last = mod(prob.N,len);
                obj.cnt(1:last) = obj.cnt(1:last) + 1;
                for i=1:len
                    if startsWith(OP{i}{1}, 'OP_SDM')
                        obj.op{i}{2} = {prob.lower,prob.upper,OP{i}{2}{1}};
                    elseif strcmp(OP{i}{1},'OP_PM')
                        obj.op{i}{2} = {prob.lower,prob.upper,OP{i}{2}{1},OP{i}{2}{2}};
                    end
                end
            end
        end
        function Offspring = Generate(obj,Pop,Prob)
        %Generate - Generate offspring from Pop and limit them to the Prob bounds.
        %
        %   Example:
        %       HC.Generate(Population,Problem)
            decs = Pop.decs;
            chaos = randperm(Prob.N);
            len = size(obj.op,2);
            if len ~= 0
                tmp = decs(chaos,:);
                s = 1;
                for i=1:len
                    e = s + obj.cnt(i) - 1;
                    decs(s:e,:) = feval(obj.op{i}{1},tmp(s:e,:),obj.op{i}{2});
                    s = e + 1;
                end
            end
            decs = max(decs,repmat(Prob.lower,Prob.N,1));
            decs = min(decs,repmat(Prob.upper,Prob.N,1));
            Offspring = Prob.Evaluation(decs);
        end
    end
end
