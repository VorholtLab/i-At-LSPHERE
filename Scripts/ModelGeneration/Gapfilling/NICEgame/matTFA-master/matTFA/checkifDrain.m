function model=checkifDrain(model)

for i=1:length(model.rxns)
    if length(find(model.S(:,i)))==1
        model.isDrain(i,1)=1;
    else
        model.isDrain(i,1)=0;
    end
end
model.isDrain(find(model.c))=1;