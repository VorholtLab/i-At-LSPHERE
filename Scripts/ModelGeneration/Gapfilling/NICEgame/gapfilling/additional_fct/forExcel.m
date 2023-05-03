function forExcel

forExcel = cell(length(active_rxns),1);
d = 1;
for i = 1:size(active_rxns,1)
    for j = 1:size(active_rxns{i},1)
        if j == 1
            if ~isempty(rxnsToCheck)
                forExcel{d,1} = rxnsToCheck{i};
                forExcel{d,2} = model.rxnNames(ismember(model.rxns,rxnsToCheck{i}));
                bbbs{1} = gfres{i,2};
            else
                forExcel{d,1} = rxnsToCheck{i};
                forExcel{d,2} = model.metNames(ismember(model.mets,rxnsToCheck{i}));
                bbbs{1} = gfres{i,2};
            end
            bbbs_name = '';
            bbbs_='';
            for s = 1:length(bbbs)
                if s == 1
                    bbbs_ = bbbs{s};
                    bbbsname = model.metNames(ismember(model.mets, bbbs{s}));
                    bbbs_name = bbbsname;
                else
                    bbbs_ = strcat(bbbs_,',',{' '},bbbs{s});
                    bbbsname = model.metNames(ismember(model.mets, bbbs{s}));
                    bbbs_name = strcat(bbbs_name,',',{' '},bbbsname);
                end
            end
            forExcel{d,3} = bbbs_;
            forExcel{d,4} = bbbs_name;
            forExcel{d,5} = 'ALTNB:1';
            if ~isempty(rxnsToCheck)
                forExcel{d,6} = printRxnFormula(model, rxnsToCheck{i});
                forExcel{d,8} = printRxnFormula(model, rxnsToCheck{i},false,false,true);
            end
            forExcel{d,7} = 'ALTNB:1';
            forExcel{d,9} = 'ALTNB:1';
            forExcel{d,10} = gfnew{i,4}{j,1};
            forExcel{d,11} = gfnew{i,4}{j,2};
            d = d + 1;
            for n = 1:length(active_rxns{i}{j,1})
                forExcel{d,5} = active_rxns{i}{j,1}{n};
                forExcel{d,6} = dir{i}{j}{n};
                forExcel{d,7} = active_rxns{i}{j,2}{n};
                forExcel{d,8} = dir{i}{j}{n};
                forExcel{d,9} = active_rxns{i}{j,3}{n};
                d = d + 1;
            end
        else
            forExcel{d,5} = strcat('ALTNB:',num2str(j));
            forExcel{d,7} = strcat('ALTNB:',num2str(j));
            forExcel{d,9} = strcat('ALTNB:',num2str(j));
            forExcel{d,10} = gfnew{i,4}{j,1};
            forExcel{d,11} = gfnew{i,4}{j,2};
            d = d + 1;
            for n = 1:length(active_rxns{i}{j,1})
                forExcel{d,5} = active_rxns{i}{j,1}{n};
                forExcel{d,6} = dir{i}{j}{n};
                forExcel{d,7} = active_rxns{i}{j,2}{n};
                forExcel{d,8} = dir{i}{j}{n};
                forExcel{d,9} = active_rxns{i}{j,3}{n};
                d = d + 1;
            end
        end
    end
end
if (isempty(rxnsToCheck) && not(isempty(forExcel)))
    for i=3:7
        forExcel(:,i)=forExcel(:,i+2);
    end
end

if gfbiomass && isempty(rxnsToCheck)
    header = {'rescued rxn', 'rescued rxn name',...
        'solutions', 'directionality',...
        'solutions formula','directionality','solutions with metnames','lump rxn of the solution', 'lump rxn with metnames'};
    forExcel(:,[8,9]) = [];
else
    if not(isempty(rxnsToCheck))
        header = {'rescued rxn keggid', 'rescued rxn name', 'bbb involved', 'bbb names', 'solutions',...
            'rescued rxn formula / direction of solution rxns', 'solution formula', 'rescued rxn with met names / direction of sol rxns',...
            'solution with met names', 'lump rxn of the solution', 'lump rxn with metnames'};
    else
        header = {'rescued bbb keggid', 'rescued bbb name', 'solutions',...
            'directionality', 'solutions formula', 'directionality',...
            'solution with metnames', 'lump rxn of the solutions', 'lump rxn with metnames'};
        forExcel(:,[8,9]) = [];
    end
end
if isempty(DPs)
    forExcel = header;
else
    forExcel = [header; forExcel];
end

save(name, 'active_rxns', 'DPs_all', 'gf_results', 'forExcel')