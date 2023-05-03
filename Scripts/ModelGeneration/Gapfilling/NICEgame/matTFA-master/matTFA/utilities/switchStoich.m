function model = switchStoich(model,indeces)
%simply switches the stoichiometric direction of reactions indicated by
%their index. MASS BALANCE MODEL ONLY

model.S(:,indeces) = -model.S(:,indeces);


ans = model.lb(indeces);
model.lb(indeces) = -model.ub(indeces);
model.ub(indeces) = -ans;