function [growthRatesAlt,predGrowth] = growAltSolnModel(ActRxns,modelDB,NumAltBounds,minMed,vitamins,nutrients,modelOrig,thermoBounds,thermoData)

[growthRatesAlt,predGrowth] = deal(zeros(length([NumAltBounds(1):NumAltBounds(2)]),length(nutrients)));

for j = 1:length(nutrients)

    model = modelOrig;

    % Define the model's medium
    model.lb(find(findExcRxns(model))) = 0;
    model.lb(find(ismember(model.rxns,intersect(model.rxns(find(findExcRxns(model))),findRxnsFromMets(model,[minMed;vitamins]))))) = -1000;
    model.lb(find(ismember(model.rxns,intersect(model.rxns(find(findExcRxns(model))),findRxnsFromMets(model,nutrients{j}))))) = -10;
    
    if thermoBounds
        definedMediumModel = getTFAbounds(model,thermoData);
    else
        definedMediumModel = model;
    end

    for n = NumAltBounds(1):NumAltBounds(2)

        model = definedMediumModel;

        for r = 1:length(ActRxns{1,1}{n,1})

            s = split(ActRxns{1,1}{n,2}{r},' ');
            metsToAdd = setdiff(s(find(contains(s,'_'))),model.mets);

            for m = 1:length(metsToAdd)
                metIndexModelDB = find(ismember(modelDB.mets,metsToAdd{m}));
                model.mets = [model.mets;metsToAdd{m}];
                model.metNames = [model.metNames;modelDB.metNames(metIndexModelDB)];
                model.metFormulas = [model.metFormulas;modelDB.metFormulas(metIndexModelDB)];
                model.metCharge = [model.metCharge;modelDB.metCharge(metIndexModelDB)];
                model.b = [model.b;modelDB.b(metIndexModelDB)];
                model.S = [model.S;zeros(1,length(model.rxns))];
                model.csense = [model.csense;'E'];
            end

            model = addReaction(model,['R_' ActRxns{1,1}{n,1}{r}],'reactionFormula',ActRxns{1,1}{n,2}{r});
        end

        FBAsoln = optimizeCbModel(model);

        growthRatesAlt(n,j) = FBAsoln.f;
    end
end
% predGrowth(find(growthRatesAlt)) = 1;
predGrowth(find(growthRatesAlt > 1e-3)) = 1;