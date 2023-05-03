function ACHRSampler_LCSB(model,warmupPoints,fileName,nFiles,pointsPerFile,stepsPerPoint,Nv_tol,initPoint,fileBaseNo,maxTime)
% ACHRSampler Artificial Centering Hit-and-Run sampler
%
% ACHRSampler(model,warmupPoints,fileName,nFiles,pointsPerFile,stepsPerPoint,initPoint,fileBaseNo,maxTime)
%
%INPUTS
% model    Model structure
% warmupPoints Warmup points
% fileName Base fileName for saving results
% nFiles   Number of files created
% pointsPerFile Number of points per file saved
% stepsPerPoint Number of sampler steps per point saved
%
%OPTIONAL INPUTS
% initPoint     Initial point (Default = centerPoint)
% fileBaseNo    Base file number for continuing previous sampler run 
%               (Default = 0)
% maxTime       Maximum time limit (Default = 36000)
%
% Markus Herrgard, Gregory Hannum, Ines Thiele, Nathan Price 4/14/06

warning off MATLAB:divideByZero;

% Set the minimum allowed Nv=0 Tolerance
if (nargin < 7)
    Nv_tol = 1e-8;
end
    
if (nargin < 9)
  fileBaseNo = 0;
else
    if (isempty(fileBaseNo))
        fileBaseNo = 0;
    end
end
if (nargin < 10)
    maxTime = 10*3600;
end

%N = null(full(model.S));

% Minimum allowed distance to the closest constraint
maxMinTol = 1e-9;
% Ignore directions where u is really small
uTol = 1e-9; 
% Project out of directions that are too close to the boundary
dTol = 1e-14;

% Number of warmup points
[nRxns,nWrmup] = size(warmupPoints);

% Find the center of the space
centerPoint = mean(warmupPoints,2);

% Set the start point
if (nargin < 8)
    prevPoint = centerPoint;
else
    if (~isempty(initPoint))
        prevPoint = initPoint;
    else
        prevPoint = centerPoint;
    end
end

fidErr = fopen('ACHRerror.txt','w');

totalStepCount = 0;
totalAcceptStepCount = 0;

h = waitbar(0,'ACHR sampling in progress ...');
totalCount = nFiles*pointsPerFile*stepsPerPoint;

t0 = cputime;
fprintf('File #\tPoint #\tStep #\tTime\t#Time left\n');
for i = 1:nFiles

    % Allocate memory for all points
    points = zeros(nRxns,pointsPerFile); 
    
    pointCount = 1;
    while (pointCount <= pointsPerFile)
            
        % Create the random step size vector
        randVector = rand(stepsPerPoint,1);
        
        stepCount = 1;
        while (stepCount <= stepsPerPoint)
            
            % Pick a random warmup point
            randPointID = ceil(nWrmup*rand);
            randPoint = warmupPoints(:,randPointID);
            
            % Get a direction from the center point to the warmup point
            u = (randPoint-centerPoint);
            u = u/norm(u);
            
            % Figure out the distances to upper and lower bounds
            distUb = (model.ub - prevPoint);
            distLb = (prevPoint - model.lb);
            
            % Figure out if we are too close to a boundary
            validDir = ((distUb > dTol) & (distLb > dTol));
            %model.rxns(~validDir)
            
            % Zero out the directions that would bring us too close to the
            % boundary. This may cause problems.
            %u(~validDir) = 0;
            
            % Figure out positive and negative directions
            posDirn = find(u(validDir) > uTol);
            negDirn = find(u(validDir) < -uTol);
            
            % Figure out all the possible maximum and minimum step sizes
            maxStepTemp = distUb(validDir)./u(validDir);
            minStepTemp = -distLb(validDir)./u(validDir);
            maxStepVec = [maxStepTemp(posDirn);minStepTemp(negDirn)];
            minStepVec = [minStepTemp(posDirn);maxStepTemp(negDirn)];
            
            % Figure out the true max & min step sizes
            maxStep = min(maxStepVec);
            minStep = max(minStepVec);
            %fprintf('%f\t%f\n',minStep,maxStep);
            
            % Find new direction if we're getting too close to a constraint
            if (abs(minStep) < maxMinTol & abs(maxStep) < maxMinTol) | (minStep > maxStep)
                fprintf('Warning %f %f\n',minStep,maxStep);
                continue;
            end
            
            % Pick a rand out of list_of_rands and use it to get a random
            % step distance
            stepDist = randVector(stepCount)*(maxStep-minStep)+minStep;
            
            % Advance to the next point
            curPoint_tmp = prevPoint + stepDist*u;
            deviation_ACHR = max(abs(model.S*curPoint_tmp-model.b));
            
            
            % Print out errors
            if (mod(totalStepCount,2000)==0)
              fprintf(fidErr,'%10.8f\t%10.8f\t',max(curPoint_tmp-model.ub),max(model.lb-curPoint_tmp));
            end

            timeElapsed = cputime-t0;
            
            % Print step information
            if (mod(totalStepCount,5000)==0)  
              timePerStep = timeElapsed/totalStepCount;
              fprintf('%d\t%d\t%d\t%8.2f\t%8.2f\n',i,pointCount,totalStepCount,timeElapsed/60,(totalCount-totalStepCount)*timePerStep/60);
            end
            
            overInd = find(curPoint_tmp > model.ub);
            underInd = find(curPoint_tmp < model.lb);
            
            if (any((model.ub-curPoint_tmp) < 0) || any((curPoint_tmp-model.lb) < 0))
              curPoint_tmp(overInd) = model.ub(overInd);
              curPoint_tmp(underInd) = model.lb(underInd);
            
            end
            
            if (mod(totalStepCount,2000)==0)
              fprintf(fidErr,'%10.8f\n',full(max(max(abs(model.S*curPoint_tmp-model.b)))));
            end
            
            % Check if curPoint_tmp satisfied Nv_tol
            if (mod(stepCount,10)==0) || (mod(stepCount,stepsPerPoint)==0)
                if (deviation_ACHR < Nv_tol)
                    curPoint = curPoint_tmp;
                    prevPoint = curPoint;
                    stepCount = stepCount + 1;
                else                    
                    randPointID = ceil((nWrmup-1)*rand);
                    prevPoint   = warmupPoints(:,randPointID);
                end
            else                
                prevPoint = curPoint_tmp;
                stepCount = stepCount + 1;
            end
            
            
            % Count the total number of steps
            totalStepCount = totalStepCount + 1;
            if mod(totalStepCount,50) == 0
                waitbar(totalStepCount/totalCount,h);
            end
            
            %recalculate the center point ONLY if the current point meets Nv_tol
            if (deviation_ACHR < Nv_tol)
                totalAcceptStepCount = totalAcceptStepCount + 1;            
                centerPoint = ((nWrmup+totalAcceptStepCount)*centerPoint + curPoint_tmp)/(nWrmup+totalAcceptStepCount+1);
            end
            
            % Exit if time exceeded
            if (timeElapsed > maxTime)
                points(:,pointCount) = curPoint;
                file = [fileName '_' num2str(fileBaseNo+i) '.mat'];
                save (file,'points');
%                 save ACHR_last_point.mat curPoint
                close(h);
                return;
            end
            
        end % Steps per point
        
        % Add the current point to points
        points(:,pointCount) = curPoint;
       
        pointCount = pointCount + 1;
         
    end % Points per cycle
    
    % Save current points to a file
    file = [fileName '_' num2str(fileBaseNo+i) '.mat'];
    save (file,'points');
    
end
close(h);

% Save last point for restart
% save ACHR_last_point.mat curPoint

fclose(fidErr);