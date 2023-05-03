function warmupPts= createHRWarmup_LCSB(model,nPoints,Nv_tol,verbFlag,nPointsCheck)
% createHRWarmup Create a warmup point set for hit-and-run sampling by
% combining orthogonal and random points
%
% warmupPts= createHRWarmup(model,nPoints,verbFlag)
%
%INPUTS
% model     Model structure
%
%OPTIONAL INPUTS
% nPoints   Number of warmup points (Default = 5000);
% verbFlag  Verbose flag (Default = false)
%
%OUTPUT
% warmupPts Set of warmup points
%
% Markus Herrgard 4/21/06
%
% Richard Que (11/23/09) Integrated subfunctions into script.

if (nargin < 2)||isempty(nPoints), nPoints = 5000; end
if (nargin < 3)||isempty(Nv_tol), Nv_tol = 1e-8; end
if (nargin < 4)||isempty(verbFlag), verbFlag = false; end
if (nargin < 5)||isempty(nPointsCheck), nPointsCheck = true; end

if isfield(model,'S')
    [nMets,nRxns] = size(model.S);
    model.A=model.S;
elseif isfield(model,'A')
    [nMets,nRxns] = size(model.A);
end
if ~isfield(model,'csense')
    model.csense(1:size(model.S,1)) = 'E';
end

if nPointsCheck && (nPoints < nRxns*2) 
    warning(['Need a minimum of ' num2str(nRxns*2) ' warmup points']);
    nPoints = nRxns*2;
end
warmupPts = sparse(nRxns,nPoints);

i = 1;
i_store = 0;
h = waitbar(0,'Creating warmup points ...');
%Generate the points
while i_store < nPoints
    if mod(i,10) == 0
        waitbar(i_store/nPoints,h);
    end
    
    % Create random objective function
    model.c = rand(nRxns,1)-0.5;
    
    for maxMin = [1, -1]
        % Set the objective function
        if i <= nRxns
            model.c = zeros(nRxns,1);
            model.c(i) = 1;
        end
        model.osense = maxMin;
        
        % Determine the max or min for the rxn
        sol = solveCobraLP(model);
        x = sol.full;
        status = sol.stat;
        if status == 1
            validFlag = true;
        else
            display ('invalid solution')
            validFlag = false;
            display(status)
            pause;
        end
        
        % Continue if optimal solution is found
        
        % Move points to within bounds
        x(x > model.ub) = model.ub(x > model.ub);
        x(x < model.lb) = model.lb(x < model.lb);
        
        deviation_wmpPts = max(abs(model.S*x-model.b));
        
        % Store point only if N*v Tolerance is satisfied
        if (maxMin == 1) 
            if (deviation_wmpPts < Nv_tol) && validFlag
                i_store = i_store+1;
                warmupPts(:,i_store) = x;
            else
                display ('N*v Tolerance not satisfied')
                fprintf('%f\n',i_store/i/2);
            end    
        else
            if (deviation_wmpPts < Nv_tol) && validFlag && i_store < nPoints   
                i_store = i_store+1;
                warmupPts(:,i_store) = x;            
            else
                display ('N*v Tolerance not satisfied')
                fprintf('%f\n',i_store/i/2);
            end
        end
        
        if (verbFlag)
            if mod(i,100)==0
                fprintf('%4.1f\n',i_store/nPoints*100);
            end
        end
        
        
    end
    if validFlag
        i = i+1;
    end 
end
centerPoint = mean(warmupPts,2);
% Move points in
warmupPts = warmupPts*.33 + .67*centerPoint*ones(1,nPoints);

close(h);
