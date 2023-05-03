function isDrain=checkDrainRxns(model)

if isfield(model,'isDrain')
    model = rmfield(model,'isDrain');
end

isDrain=zeros(length(model.rxns),1);

for i=1:length(model.rxns)
    if nnz(model.S(:,i)) == 1
        isDrain(i,1)=1;
    end
end
