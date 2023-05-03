function [active_rxns, gf_results, DPs_all, forExcel] = gapFillingBBB(model,NumAlt,tagGFbiom,rxnsToCheck,bbbsToCheck,tagMin,time,name)

%% Description
%   perform a gapfilling analysis on a merged model from PrepareForGapFilling
%
% INPUTS
%   model                   a merged model with the function PrepareForGapFilling
%   NumAlt                  number of alternative solutions to find per rxn
%                           and/or bbb, opt default 1
%   tagGFbiom               gap fill for biomass and not for bbbs, opt
%                           default true
%   rxnsToCheck             rxns to KO and gap fill (superessentiality
%                           studies), opt default all rxns in model.rxnRescued
%   bbbsToCheck             BBBs to gap fill for biomass or for each rxn in
%                           rxnsToCheck, opt default bbbs in model.bbbNotProduced
%   name                    name of the file in which the data will be saved,
%                           the script will save the data at each iteration


% OUTPUTS
%   active_rxns             contains the rxns of each gapfilling solutions
%   overall_am_names_full   contains the bbbs not produced for each rxns tested
%   DPs_all                 contains the solution vector of each gapfilling
%   rxn_can_produce         contains the rxns (or bbb) for which the model didn't have to incorporate any rxns from the DB model
%   rxn_cannot_produce      contains the rxns for which the model didn't find a solution  a solution
%   gf_results              results of the gapfilling ordered so it is easy to read and use:
%                           column 1: rxn rescued
%                           column 2: rxn that fill the gap
%                           column 3: bbbs saved by gapfilling
%
%   forExcel                contains the results ordered for Excel, just copy/paste it to excel
%
%   Yannick Francioli & Anush Chiappino-Pepe 2017

if (nargin < 2)
    NumAlt = 1;
end
if (nargin < 3)
    tagGFbiom = 1;
end
if (nargin < 4)
    rxnsToCheck = model.rxnRescued;
else
    if not(isempty(rxnsToCheck))
        rxnsToCheck = model.rxnRescued(ismember(model.rxnRescued,rxnsToCheck));
    end
end
if (nargin < 5)
    bbbsToCheck = model.bbbRescued;
% else
%     if not(isempty(rxnsToCheck))
%         rxnsToCheck = model.rxnRescued(ismember(model.rxnRescued,rxnsToCheck));
%     end
end
if (nargin < 6)
    tagMin = 1;
end
if (nargin < 7)
    time = 2000;
end
if (nargin < 8)
    name = strcat('gapFilling','_',model.id);
end

biomassRxn = model.rxns(model.c==1);
if isempty(rxnsToCheck)
    tagGFbiom = true;
end
if tagGFbiom
    growthReq = 0.08;
    fprintf(strcat('gap-filling per biomass with growth requirement of ',num2str(growthReq),'\n'));
else
    growthReq = 0.1*model.growthWT;
    fprintf(strcat('gap-filling per BBBs with growth requirement of ',num2str(growthReq),'\n'));
end

% check that there is no objective function in the database model
model.c(model.rxnIndDB==1)=0;

% no possible combinations
if ~tagGFbiom && isempty(rxnsToCheck)
    error('provide true tagGFbiom or a set of rxnsToCheck')
end
if ~tagGFbiom && isempty(bbbsToCheck)
    error('provide true tagGFbiom or a set of bbbsToCheck')
end

%% Performing GF
% block all demand reactions and the biomass rxn
fprintf('performing gap-filling\n');
model.var_ub(ismember(model.varNames,strcat('F_DM_',strrep(model.mets(model.S(:,model.c==1)<0),'-','_'))))=0;
model.var_ub(ismember(model.varNames,strcat('R_DM_',strrep(model.mets(model.S(:,model.c==1)<0),'-','_'))))=0;
model.var_ub(ismember(model.varNames,strcat('F_', biomassRxn)))=0;
bbb_not={'h2o_c', 'C00001_c'};
SSSS=full(model.S);
DPs=[];
arxns={};
if ~isempty(rxnsToCheck)
    % Preallocate vectors
    gfres=cell(length(rxnsToCheck),1);
    active_rxns=cell(length(rxnsToCheck),1);
    DPs_all=cell(length(rxnsToCheck),1);
    
    for i=1:length(rxnsToCheck)
        fprintf(strcat('performing gap-filling for rxn ',rxnsToCheck{i},'\n'));
        GFmodel = model;
        %knockout the essential rxn
        GFmodel.var_ub(ismember(GFmodel.varNames, strcat( 'F_',rxnsToCheck{i})))=0;
        GFmodel.var_ub(ismember(GFmodel.varNames, strcat( 'R_',rxnsToCheck{i})))=0;
        if tagGFbiom
            GFmodel.var_ub(ismember(GFmodel.varNames,strcat('F_', biomassRxn))) = 50;
            GFmodel.var_lb(ismember(GFmodel.varNames,strcat('F_', biomassRxn))) = growthReq;
            sol=optimizeThermoModel(GFmodel);
            
            if ~isempty(sol.x)
                gfres{i,2}='biomass';
                rxnsBFUSE=find(sol.x(ind_bfuse)==0,1); % find active use variables % rows for the rxns added from model2 (BFUSE=0)
                if isempty(rxnsBFUSE) % if no rxns were added from model2 the model to gap-fill could produce this specific BBB
                    gfres{i,1}='can_be_produced';
                else
                    gfres{i,1}='gap_filled';
                    fprintf(strcat('getting alternative solutions\n'));
                    [DPs,~,~] = findDP4gf(GFmodel,NumAlt,indUSE,time,tagMin); %change Num allternative
                    if ~isempty(DPs) % loop to find min+1
                        GFmodel.rhs(ismember(GFmodel.constraintNames,'CUT_0'))=0;
                        [DPs2,~,~] = findDP4gf(GFmodel,NumAlt,indUSE,time,tagMin);
                        if ~isempty(DPs2)
                            DPs = [DPs, DPs2];
                        end
                    end
                end
            else
                gfres{i,1}='cannot be gap-filled';
            end
            
            if ~isempty(DPs)
                arxns=cell(size(DPs,2),3);
                for j=1:size(DPs,2)
                    bfuse=DPs(ind_bfuse,j);
                    rxnfor=strrep(GFmodel.varNames(ind_bfuse(bfuse==0)),'BFUSE_', '');
                    rxnfor=strrep(rxnfor, 'BU_', '');
                    [~, bpath2]=ismember(rxnfor, GFmodel.rxns);
                    if ~isempty(bpath2)
                        arxns{j,1}=GFmodel.rxns(bpath2);
                        arxns{j,2}=printRxnFormula(GFmodel, GFmodel.rxns(bpath2));
                        arxns{j,3}=printRxnFormula(GFmodel, GFmodel.rxns(bpath2),false,false,true);
                    end
                end
            end
            active_rxns{i} = arxns;
            DPs_all{i} = DPs;
            
        else
            for bbbj=1:length(bbbsToCheck{i})
                if ~ismember(bbbsToCheck{i}, bbb_not) % for all BBBs (not h2o)
                    GFmodeli=GFmodel;
                    % extract info about the stoichiometric coefficients in the biomass rxn
                    stoich_bbb=SSSS(ismember(GFmodel.mets,bbbsToCheck{i}{bbbj}), GFmodel.c==1);
                    GFmodeli.var_ub(ismember(GFmodel.varNames,strcat('F_DM_', strrep(bbbsToCheck{i}{bbbj}, '-', '_'))))=50; % do not allow uptake of BBB
                    GFmodeli.var_lb(ismember(GFmodel.varNames,strcat('F_DM_', strrep(bbbsToCheck{i}{bbbj}, '-', '_'))))=-stoich_bbb*growthReq;
                    sol=optimizeThermoModel(GFmodeli);
                else
                    sol.x=[];
                end
                % store solutions per BBB
                gfres{i,2}{bbbj}=bbbsToCheck{i}{bbbj};
                if ~isempty(sol.x)
                    rxnsBFUSE=find(sol.x(ind_bfuse)==0,1); % find active use variables % rows for the rxns added from model2 (BFUSE=0)
                    if isempty(rxnsBFUSE) % if no rxns were added from model2 the model to gap-fill could produce this specific BBB
                        gfres{i,1}{bbbj}='can_be_produced';
                    else
                        gfres{i,1}{bbbj}='gap_filled';
                        fprintf(strcat('getting alternative solutions\n'));
                        [DPs,~,~] = findDP4gf(GFmodel,NumAlt,indUSE,time,tagMin); %change Num allternative
                        if ~isempty(DPs) % loop to find min+1
                            GFmodel.rhs(ismember(GFmodel.constraintNames,'CUT_0'))=0;
                            [DPs2,~,~] = findDP4gf(GFmodel,NumAlt,indUSE,time,tagMin);
                            if ~isempty(DPs2)
                                DPs = [DPs, DPs2];
                            end
                        end
                    end
                else
                    gfres{i,1}{bbbj}='cannot be gap-filled';
                    DPs=[];
                    arxns={};
                end
                if ~isempty(DPs)
                   arxns=cell(size(DPs,2),3);
                    for j=1:size(DPs,2)
                        bfuse=DPs(ind_bfuse,j);
                        rxnfor=strrep(GFmodel.varNames(ind_bfuse(bfuse==0)),'BFUSE_', '');
                        rxnfor=strrep(rxnfor, 'BU_', '');
                        [~, bpath2]=ismember(rxnfor, GFmodel.rxns);
                        if ~isempty(bpath2)
                            arxns{j,1}=GFmodel.rxns(bpath2);
                            arxns{j,2}=printRxnFormula(GFmodel, GFmodel.rxns(bpath2));
                            arxns{j,3}=printRxnFormula(GFmodel, GFmodel.rxns(bpath2),false,false,true);
                        end
                    end
                end
                %active_rxns{i}{bbbj}=arxns;
                if isempty(arxns)
                    active_rxns{i}{bbbj}={};
                else
                active_rxns{i}{bbbj,1}=arxns{:,1};
                active_rxns{i}{bbbj,2}=arxns{:,2};
                active_rxns{i}{bbbj,3}=arxns{:,3};
                end
                DPs_all{i}{bbbj}=DPs;
            end
        end 
    end
else
    if ~tagGFbiom
        error('for no gf per biomass a set of rescued rxns should be provided')
    else
        fprintf('performing gap-filling for biomass\n');
        gfres=cell(1,2);
        active_rxns=cell(1,1);
        DPs_all=cell(1,1);
        GFmodel = model;
        GFmodel.var_ub(ismember(GFmodel.varNames,strcat('F_', biomassRxn))) = 50;
        GFmodel.var_lb(ismember(GFmodel.varNames,strcat('F_', biomassRxn))) = growthReq;
        sol=optimizeThermoModel(GFmodel);
    end
    % store solution for biomass
    if ~isempty(sol.x)
        i=1;
        gfres{i,2}='biomass';
        rxnsBFUSE=find(sol.x(ind_bfuse)==0,1); % find active use variables % rows for the rxns added from model2 (BFUSE=0)
        if isempty(rxnsBFUSE) % if no rxns were added from model2 the model to gap-fill could produce this specific BBB
            gfres{i,1}='can_be_produced';
        else
            gfres{i,1}='gap_filled';
            fprintf(strcat('getting alternative solutions\n'));
            [DPs,~,~] = findDP4gf(GFmodel,NumAlt,indUSE,time,tagMin); %change Num allternative
            if ~isempty(DPs) % loop to find min+1
                GFmodel.rhs(ismember(GFmodel.constraintNames,'CUT_0'))=0;
                [DPs2,~,~] = findDP4gf(GFmodel,NumAlt,indUSE,time,tagMin);
                if ~isempty(DPs2)
                    DPs = [DPs, DPs2];
                end
            end
        end
    else
        gfres{i,1}='cannot be gap-filled';
        DPs=[];
        arxns=[];
    end    
    if ~isempty(DPs)
        arxns=cell(size(DPs,2),3);
        for j=1:size(DPs,2)
            bfuse=DPs(ind_bfuse,j);
            rxnfor=strrep(GFmodel.varNames(ind_bfuse(bfuse==0)),'BFUSE_', '');
            rxnfor=strrep(rxnfor, 'BU_', '');
            [~, bpath2]=ismember(rxnfor, GFmodel.rxns);
            if ~isempty(bpath2)
                arxns{j,1}=GFmodel.rxns(bpath2);
                arxns{j,2}=printRxnFormula(GFmodel, GFmodel.rxns(bpath2));
                arxns{j,3}=printRxnFormula(GFmodel, GFmodel.rxns(bpath2),false,false,true);
            end
        end
    end
    active_rxns{i} = arxns;
    DPs_all{i} = DPs; 
end

gf_results = cell(length(active_rxns),1);
for k = 1:length(active_rxns)
    if tagGFbiom && isempty(rxnsToCheck)
        gf_results{k,1} = 'biomass';
        gf_results{k,3} = gfres{k,1};
    else
        gf_results{k,1} = rxnsToCheck{i};
        gf_results{k,3} = gfres{i,1};
    end
    for n = 1:size(active_rxns{k},1)
       
            gf_results{k,2}{n,1} = active_rxns{k}{n,1};
     
        
            gf_results{k,2}{n,2} = active_rxns{k}{n,2};
 
      
            gf_results{k,2}{n,3} = active_rxns{k}{n,3};
      
    end
end
save(name, 'active_rxns', 'DPs_all', 'gf_results')
fprintf('gap-filling done\n');
%% get directionality informations
fprintf('getting directionality information for gap-filling results\n');
dir = {};%cell(length(active_rxns),1);
DP={};
for i = 1:length(active_rxns)
    for j = 1:size(active_rxns{i},1)
        for n = 1:length(active_rxns{i}{j,1})
            f = strcat('F_',active_rxns{i}{j,1}{n});
            b = strcat('R_',active_rxns{i}{j,1}{n});
            DP = DPs_all{i}{:,j};
            if not(isempty(DP))
            if DP((ismember(model.varNames,f))) > 0
                 gf_results{i,2}{j,4}{n}= 'forward';
                dir{i}{j}{n} = 'forward';
            end
            if DP((ismember(model.varNames,b))) > 0
                gf_results{i,2}{j,4}{n} = 'backward';
                dir{i}{j}{n} = 'backward';
            end
            else
                 gf_results{i,2}{j,4}{n}= '';
                 dir={};
            end
        end
    end
end

%% get the lump rxn for each alternative
fprintf('getting lumped rxns for gap-filling results\n');
gfnew = gf_results;
for i = 1:size(gf_results,1)
    for j = 1:size(gf_results{i,2},1)
        nn = length(gf_results{i,2}{j,1});
        if nn == 1;
            gfnew{i,4}{j,1} = printRxnFormula(model,gf_results{i,2}{j,1}{nn});
            gfnew{i,4}{j,2} = printRxnFormula(model,gf_results{i,2}{j,1}{nn},false,false,true);
        else
            all_s = zeros(size(model.S,1),nn);
            for n = 1:nn
                s = model.S(:,ismember(model.rxns,gf_results{i,2}{j,1}{n}));
                if strcmp(gf_results{i,2}{j,4}{n},'backward')
                    s = s*-1;
                end
                all_s(:,n) = s;
            end
            s_lump = sum(all_s,2);
            %s_lump(find(s_lump));
            rct_ind = find(s_lump < 0);
            prod_ind = find(s_lump > 0);
            rctmets = model.mets(rct_ind);
            rctnames = model.metNames(rct_ind);
            prodmets = model.mets(prod_ind);
            prodnames = model.metNames(prod_ind);
            rct = '';
            rctn = '';
            prod = '';
            prodn = '';
            for k = 1:length(rct_ind)
                if k == 1
                    rct = rctmets{k};
                    rctn = rctnames{k};
                else
                    rct = strcat(rct,{' + '},rctmets{k});
                    rctn = strcat(rctn,{' + '},rctnames{k});
                end
            end
            for k = 1:length(prod_ind)
                if k == 1
                    prod = prodmets{k};
                    prodn = prodnames{k};
                else
                    prod = strcat(prod,{' + '},prodmets{k});
                    prodn = strcat(prodn,{' + '},prodnames{k});
                end
            end
            gfnew{i,4}{j,1} = strcat(rct, {' <=> '}, prod);
            gfnew{i,4}{j,2} = strcat(rctn, {' <=> '}, prodn);
        end
    end
end
gf_results = gfnew;
%% put the results in a away that easy to read and that can be pasted in excel
fprintf('storing gap-filling results in an excel friendly format\n');
if tagGFbiom && isempty(rxnsToCheck)
    rxnsToCheck = biomassRxn;
end

end

