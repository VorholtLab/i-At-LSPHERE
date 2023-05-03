function [DPs,model] = findAltBBB(model,numAlt,indUSE)

% Description
% Function called inside "findBBBforStoich.m"
% This function identifies alternative solutions for the MILP problem
% defined in "findBBBforStoich.m"

if (nargin < 2)
    numAlt = 1;
end
if (nargin < 3)
    indUSE = getAllVar(model,{'BBB'});
end

NumSols = 0;
model.objtype = 1; %minimize:1, maximize:-1 (minimize for old)

sol = optimizeThermoModel(model);
% find alternate solutions
fprintf('getting alternative solutions\n');
while ((NumSols < numAlt) && ~(isempty(sol.x)))
    
    [numCons,numVars] = size(model.A);
    
    if ~(isempty(sol.x))
        NumSols = NumSols + 1;
        DPs(:,NumSols) = sol.x;
        
        USEvecSol = zeros(numVars,1);
        USEvecSol(indUSE) = sol.x(indUSE);
        actUSEvec = find(USEvecSol>0.98);
        
        NewCons = zeros(1,numVars);
        NewCons(actUSEvec) = 1;
        
        model.A(numCons+1,:) = NewCons;
        model.rhs(numCons+1) = length(actUSEvec)-0.5;
        model.constraintNames{numCons+1} = ['CUT_' num2str(NumSols)];
        model.constraintType{numCons+1} = '<'; 
        sol = optimizeThermoModel(model);
        
        if isempty(sol.x)
            break
        end
        
        fprintf('Number of DPs:\t%d\n',NumSols);
        %             if rem(NumSols,10) == 0
        %                save DPs DPs;
        %             end
    end
end

