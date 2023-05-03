function [models,conflict] = evalMetsModels4gf(sourceModel,models)

conflict = struct;
conflict.dupMet = {};

mets = sourceModel.mets;
for tt = 1:length(sourceModel.CompartmentData.compSymbolList)
    mets = strrep(mets, strcat('_', sourceModel.CompartmentData.compSymbolList{tt}), '');
end

for n = 1:numel(models)
    % Preprocess the mets from the sourcemodel to compare
    mets2 = models{n}.mets;
    for tt = 1:length(sourceModel.CompartmentData.compSymbolList)
        mets2 = strrep(mets2, strcat('_', sourceModel.CompartmentData.compSymbolList{tt}), '');
    end
    
    % Check that there are no conflicts with the metSEEDID and save the
    % conflicts
    % Correct this conflict by defining the metSEEDID from the sourcemodel
    conflict.seedid = [];
    present = 0;
    for i = 1:length(mets)
        if ismember(mets(i),mets2)
            present = present+1;
            [~,rrr] = ismember(mets(i),mets2);
            if ~isequal(sourceModel.metSEEDID(i),models{n}.metSEEDID(rrr(1)))
                conflict.seedid = [conflict.seedid; ...
                    [mets(i), sourceModel.metSEEDID(i), models{n}.metSEEDID(rrr(1))]];
                models{n}.metSEEDID(ismember(mets2,mets(i))) = sourceModel.metSEEDID(i);
            end
        end
    end
    
    % Check that there are no conflicts with the metNames and save the
    % conflicts
    % Correct this conflict by defining the metName from the sourcemodel
    conflict.metNames = [];
    present1 = 0;
    for i = 1:length(mets)
        if ismember(mets(i),mets2)
            present1 = present1+1;
            [~,rrr] = ismember(mets(i),mets2);
            if ~isequal(sourceModel.metNames(i),models{n}.metNames(rrr(1)))
                conflict.metNames = [conflict.metNames; ...
                    [mets(i), sourceModel.metNames(i), models{n}.metNames(rrr(1))]];
                models{n}.metNames(ismember(mets2,mets(i))) = sourceModel.metNames(i);
            end
        end
    end
    
    %check for Met id conflict (same name, different id)
    conflict.confMetId = {};
    k = 1;
    for i=1:length(sourceModel.metNames)
        if ismember(sourceModel.metNames(i),models{n}.metNames)
            if ~isequal(sourceModel.mets(i),models{n}.mets(ismember(models{n}.metNames,sourceModel.metNames(i))))...
                    && isequal(sourceModel.metCompSymbol(i), models{n}.metCompSymbol(ismember(models{n}.metNames,sourceModel.metNames(i))))
                
                conflict.confMetId(k,1) = sourceModel.metNames(i);
                conflict.confMetId(k,2) = sourceModel.mets(i);
                for j = 1:length(models{n}.mets(ismember(models{n}.metNames,sourceModel.metNames(i))))
                    mNames = models{n}.mets(ismember(models{n}.metNames,sourceModel.metNames(i)));
                    conflict.confMetId(k,2 + j) = mNames(j);
                end
                k = k + 1;
            end
        end
    end
end