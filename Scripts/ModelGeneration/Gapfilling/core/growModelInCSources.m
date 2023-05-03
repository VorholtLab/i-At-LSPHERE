function pG = growModelInCSources(mdl,mm,vits,nts,vmax)

    [growthRatesAlt,pG] = deal(zeros(1,length(nts)));
    for j = 1:length(nts)

        modelTest = mdl;

        % Define the model's medium
        modelTest.lb(find(findExcRxns(modelTest))) = 0;
        modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,[mm;vits]))))) = -1000;
        modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,nts{j}))))) = -vmax;

        FBAsoln = optimizeCbModel(modelTest);
        growthRatesAlt(j) = FBAsoln.f;
    end
    pG(find(growthRatesAlt > 1e-3)) = 1;
end