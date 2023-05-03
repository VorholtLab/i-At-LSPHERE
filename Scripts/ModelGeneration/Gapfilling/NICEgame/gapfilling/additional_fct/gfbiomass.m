function [resultStat,ActRxns,DPsAll] = gfbiomass(model,indUSE,NumAlt,time,tagMin,tagMinP1,tagSave,filename,DPsi)

if (nargin < 2)
    indUSE = getAllVar(model,{'BFUSE'});
    if isempty(indUSE)
        error('BFUSE integers not found')
    end
end
if (nargin < 3)
    NumAlt = 1;
end
if (nargin < 4)
    time = [];
end
if (nargin < 5)
    tagMin = 0;
end
if (nargin < 6)
    tagMinP1 = 0;
end
if (nargin < 7)
    tagSave = 0;
end
if (nargin < 8)
    filename = 'gf';
end
if (nargin < 9)
    DPsi = [];
end

resultStat = cell(1,2);
ActRxns = cell(1,1);
DPsAll = cell(1,1);
DPs = [];
arxns = {};

sol = solveTFAmodelCplex(model,300);

resultStat{1,2} = 'biomass'; % store status for biomass
if ~isempty(sol.x)
    % find active use variables: rows for the rxns added from DBmodel (BFUSE=0)
    if ~any(sol.x(indUSE)<0.1) % if no rxns were added from DBmodel the model to gap-fill could produce biomass
        resultStat{1,1} = 'can_be_produced';
    else
        resultStat{1,1} = 'gap_filled';
%         fprintf(strcat('getting alternative solutions\n'));
        [DPs,model] = findDP4gf(model,NumAlt,indUSE,time,tagMin,tagSave,filename);
        if ~isempty(DPsi) && isempty(DPs)
            DPs = DPsi;
        elseif ~isempty(DPsi) && ~isempty(DPs)
            DPs = [DPsi, DPs];
        end
        if tagSave
            save(strcat(filename,'_DPsAll.mat'), 'DPs');
        end
        if ~isempty(DPs) && tagMinP1 % loop to find min+1
            checkDif = abs(max(round(sum(DPs(indUSE,:)),0)) - min(round(sum(DPs(indUSE,:)),0)));
            if checkDif < 0.1
                model.rhs(ismember(model.constraintNames,'CUT_0')) = max(round(sum(DPs(indUSE,:)),0)) - 5.5; %limit min+1 to be 5 step different
                [DPs2,model] = findDP4gf(model,NumAlt,indUSE,time,tagMin,tagSave,filename);
                if ~isempty(DPs2)
                    DPs = [DPs, DPs2];
                end
            end
        end
        if tagSave
            save(strcat(filename,'_DPsAll.mat'), 'DPs');
        end
    end
else
    resultStat{1,1} = 'cannot be gap-filled';
end

if ~isempty(DPs)
    arxns = getActRxns(model,DPs,indUSE);
end
ActRxns{1} = arxns;
DPsAll{1} = DPs;
