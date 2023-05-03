function model = addNewConstraintInTFA(model, CName, CType, CLHS, CRHS)
% This function adds one new constraint to the model, and specifies how the
% associated variables are involved in this constraint
[num_constr,~] = size(model.A);
model.constraintNames{num_constr+1,1} = CName; % strcat('DFSEU_',model.mets{i});
model.constraintType{num_constr+1,1}  = CType; % type of constraint: '<', or '>', or '='
model.rhs(num_constr+1) = CRHS; % value of the right hand side
model.A(num_constr+1,CLHS.varIDs) = CLHS.varCoeffs; % coefficients of the involved variables
end