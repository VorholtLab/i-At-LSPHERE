function actRxns = getActRxns(model,DPs,indUSE)

actRxns = cell(size(DPs,2),3);
for j = 1:size(DPs,2)
    [y, row] = ismember(strrep(model.varNames(indUSE(DPs(indUSE,j)<0.1)),'BFUSE_', ''), model.rxns);
    if all(y)
        actRxns{j,1} = model.rxns(row);
        actRxns{j,2} = printRxnFormula(model, model.rxns(row),0);
        actRxns{j,3} = printRxnFormula(model, model.rxns(row),0,0,1);
        actRxns{j,4} = DPs(ismember(model.varNames,strrep(model.varNames(indUSE(DPs(indUSE,j)<0.1)),'BFUSE','FU')),j); %extracting directionality information
        actRxns{j,5} = getLumpRxn(model,actRxns{j,1},actRxns{j,4});
    else
        fprintf('not all rxns found!\n');
    end
end