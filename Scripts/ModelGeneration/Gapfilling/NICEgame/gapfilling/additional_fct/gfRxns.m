function [resultStat,ActRxns,DPsAll] = gfRxns(model,rxnsToCheck,indUSE,NumAlt,time,tagMin,tagSave,filename,DPsi)

if (nargin < 3)
    indUSE = getAllVar(model,{'BFUSE'});
    if isempty(indUSE)
        error('BFUSE integers not found')
    end
end
if (nargin < 4)
    NumAlt = 1;
end
if (nargin < 5)
    time = [];
end
if (nargin < 6)
    tagMin = 0;
end
if (nargin < 7)
    tagSave = 1;
end
if (nargin < 8)
    filename = 'gf';
    path = '/Users/evangelia/Documents/MATLAB/gapflling paper/gf results/'; % for cluster
end
if (nargin < 9)
    DPsi = [];
end
filename = 'gf';
    path = '/Users/evangelia/Documents/MATLAB/gapflling paper/outputs/'; 
% Preallocate vectors
resultStat = cell(length(rxnsToCheck),1);
ActRxns = cell(length(rxnsToCheck),1);
DPsAll = cell(length(rxnsToCheck),1);

for i = 1:length(rxnsToCheck)
    fprintf(strcat('performing gap-filling for rxn ',rxnsToCheck{i},'\n'));
    GFmodel = model;
    DPs = [];
    arxns = {};
    %knockout the essential rxn
    GFmodel.var_ub(ismember(GFmodel.varNames, strcat('F_',rxnsToCheck{i}))) = 0;
    GFmodel.var_ub(ismember(GFmodel.varNames, strcat('R_',rxnsToCheck{i}))) = 0;
    sol = solveTFAmodelCplex(GFmodel,300);
    
    if ~isempty(sol.x)
        resultStat{i,2} = 'biomass';
        % find active use variables % rows for the rxns added from model2 (BFUSE=0)
        if ~any(sol.x(indUSE)<0.1) % if no rxns were added from model2 the model to gap-fill, the model could produce biomass with this KO
            resultStat{i,1} = 'can_be_produced';
        else
            resultStat{i,1} = 'gap_filled';
            fprintf(strcat('getting alternative solutions for min\n'));
            [DPs,GFmodel] = findDP4gf(GFmodel,NumAlt,indUSE,time,tagMin,tagSave,filename); %change Num allternative
            if ~isempty(DPsi) && isempty(DPs)
                DPs = DPsi;
            elseif ~isempty(DPsi) && ~isempty(DPs)
                DPs = [DPsi, DPs];
            end
            if tagSave
                save(strcat(path,filename,'_DPsAll.mat'), 'DPs');
            end
            if ~isempty(DPs) % loop to find min+1
                checkDif = abs(max(round(sum(DPs(indUSE,:)),0)) - min(round(sum(DPs(indUSE,:)),0)));
                if checkDif < 0.1
                    GFmodel.rhs(ismember(GFmodel.constraintNames,'CUT_0')) = max(round(sum(DPs(indUSE,:)),0)) - 5.5;
                    fprintf(strcat('getting alternative solutions for min+1\n'));
                    [DPs2, GFmodel] = findDP4gf(GFmodel,NumAlt,indUSE,time,tagMin,tagSave,filename);
                    if ~isempty(DPs2)
                        DPs = [DPs, DPs2];
                    end
                end
            end
            if tagSave
                save(strcat(path,filename,'_DPsAll.mat'), 'DPs');
            end
        end
    else
        resultStat{i,1} = 'cannot be gap-filled';
    end
    if ~isempty(DPs)
        arxns = getActRxns(model,DPs,indUSE);
    end
    ActRxns{i} = arxns;
    DPsAll{i} = DPs;
end