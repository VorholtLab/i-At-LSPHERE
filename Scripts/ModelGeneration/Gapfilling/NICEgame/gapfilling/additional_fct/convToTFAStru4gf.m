function model = convToTFAStru4gf(model,DBThermo,flagTFA,rxnsNoThermo,metabData)

if (nargin < 2)
    DBThermo = [];
end
if (nargin < 3)
    flagTFA = 0;
end
if (nargin < 4)
    rxnsNoThermo = [];
end
if (nargin < 5)
    metabData = [];
end

if isempty(DBThermo) && flagTFA
    load DB_AlbertyUpdate_keng;
    DBThermo = DB_AlbertyUpdate;
end

sol = optimizeCbModel(model);
bbb = model.mets(model.S(:,model.c==1)<0);
model = addDemandReaction4gf(model, bbb, 0, 0, 0);

if ~flagTFA
    model = convGEMToTFAStruc(model);
else
    model = prepModelforTFA(model, DBThermo, model.CompartmentData, 0, 0);
    model = convToTFA(model, DBThermo, rxnsNoThermo,'DGo', [], 0.1*sol.f, [], [], 0);
    if ~isempty(metabData) % add metabolomics data
        model = loadConstraints(model,metabData);
    end
end
end