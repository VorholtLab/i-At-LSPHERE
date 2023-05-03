function lumpedRxn = getLumpRxn(model,rxns,FUres)

stoich = zeros(length(model.mets),length(rxns));
for i = 1:length(rxns)
    stoich(:,i) = model.S(:,ismember(model.rxns,rxns(i)));
    if FUres(i)==0 % in the TFA formulation if FU=0 the F_rxn = 0
        stoich(:,i) = stoich(:,i)*-1;
    end
end

model.S = sum(stoich,2);
lumpedRxn{1,1} = printRxnFormula(model, model.rxns(1),0);
lumpedRxn{1,2} = printRxnFormula(model, model.rxns(1),0,0,1);
