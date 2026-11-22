function Offspring = OP_SDM(Parent,Args)
% SDM operator.

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

    [lower,upper,ratio] = deal(Args{:});
    [N,D] = size(Parent);
    Offspring = Parent;
    if N ~= 0 && D ~= 0
        len = max(floor(D*ratio),1);
        sle = randi([-1,1],1,N);
        for i=1:N
            if sle(i) > 0
                ch = find(Parent(i,:) ~= 0);
                if ~isempty(ch)
                    [~,chi] = sort(sum(Parent(:,ch) ~= 0,1));
                    ch = ch(chi);
                    ch = ch(1:min(len,length(ch)));
                    Offspring(i,ch) = 0;
                end
            elseif sle(i) < 0
                ch = find(Parent(i,:) == 0);
                if ~isempty(ch)
                    [~,chi] = sort(sum(Parent(:,ch) ~= 0,1),'descend');
                    ch = ch(chi);
                    ch = ch(1:min(len,length(ch)));
                    Offspring(i,ch) = rand(1,length(ch)).*(upper(ch)-lower(ch))+lower(ch);
                end
            else
                k = randi(D);
                s = randi(D);
                if s < k
                    Offspring(i,:) = Offspring(i,[1:s-1,k,s:k-1,k+1:end]);
                elseif s > k
                    Offspring(i,:) = Offspring(i,[1:k-1,k+1:s-1,k,s:end]);
                end
            end
        end
    end
end