function Population = InitPopulation(prob, sLower, sUpper)
% Initialize the population of MOEA-EO.

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

    %% Result set up
    pop = prob.Initialization();
    varCount = size(prob.lower,2);
    mask = false(prob.N, varCount);

    %% Determine the positioning of each stripe per individual
    densityVector = 1 - linspace(sLower, sUpper, prob.N);

    widthVector = round(densityVector.*prob.D);
    
    % Put widths back into bound if rounding error occurred
    lb = floor((1- sLower)*prob.D);
    widthVector(widthVector > lb) = lb;

    cumulativeWidths = cumsum(widthVector);

    % if all sparsities are 100%, then skip processing, since everything
    % will be zeros 
    if sum(widthVector == 0) == prob.N
        processedIndvs = prob.N; 
    else
        processedIndvs = 0;
    end

    cycle_count = 0;
    cycles = zeros(prob.N, prob.D);
    while processedIndvs < prob.N
        % Figure out how many stripes will fit in this cycle 
        cycle_count = cycle_count + 1;
        spotsThatFitMask = cumulativeWidths <= prob.D & cumulativeWidths ~= 0;
        numThatFit = sum(spotsThatFitMask);
        largestFit = max(cumulativeWidths(spotsThatFitMask));

        cumulativeWidths = cumulativeWidths - largestFit; 
        cumulativeWidths(cumulativeWidths<0) = 0; 
        processedIndvs = processedIndvs + numThatFit;
        spotsThatFit = find(spotsThatFitMask);
        cycles(cycle_count,1:numThatFit) = spotsThatFit;   

    end

    %% Create density mask 

    % Mask out non-zero values cycle-by-cycle
    currentIndv = 1;
    for c = 1:cycle_count
        cycle = cycles(c, cycles(c, :) ~= 0);

        widths = widthVector(cycle);

        gapToFill = prob.D - sum(widths);
        gapSize = ceil((prob.D - sum(widths))/numel(widths));
        % For each individual in the cycle 
        position = 1;
        for i = 1:numel(widths)

            width = widths(i);

            % Determine if a gap is needed 
            gapWidth = 0;
            if gapToFill > 0
                gapWidth = gapSize;
                gapToFill = gapToFill - gapWidth;
            end

            % Determine the position of the stripe
            startPoint = position;

            if c == cycle_count
                endPoint = position+width-1; 
            else
                endPoint = position+width-1+gapWidth; 
            end

            % Prevent overflow from a gap calculation
            if endPoint > prob.D 
                endPoint = prob.D;
            end
    
            % Mask out stripe
            mask(currentIndv, startPoint:endPoint) = true; 

            % Go to the next individual 
            
            position = position+width+gapWidth;
            %position = position+width;

            currentIndv = currentIndv + 1;

        end

    end

    %% Mask off population according to stripe position 
    sparse_pop = pop.decs;
    sparse_pop(~mask) = 0;

    % Recalculate objective and constraints
    popDec = prob.CalDec(sparse_pop);
    popObj = prob.CalObj(sparse_pop);
    popCon = prob.CalCon(sparse_pop);

    Population = SOLUTION(popDec, popObj, popCon);
end