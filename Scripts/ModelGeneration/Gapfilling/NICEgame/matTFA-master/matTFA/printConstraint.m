function print_out_vec = printConstraint(model,constraintNames,ScreenPrintflag)
% OUTPUTS
% - print_out_vec: A 1 column cell of strings with the all the constraint
%   expressions
% 
% INPUTS
% 
% - model:
% - constraintNames:
% - ScreenPrintflag: Print the output on the screen or not (default true).

if nargin<3
    ScreenPrintflag = true;
end

const_indices = find(ismember(model.constraintNames,constraintNames));
A = model.A;
varNames = model.varNames;
constraintNames = model.constraintNames;
constraintType = model.constraintType;
rhs = model.rhs;

print_out_vec = [];

for j=1:length(const_indices)
    index = const_indices(j);
    variables = varNames(find(A(index,:)));
    coefficients = full(A(index,find(A(index,:))));
    constraintName = constraintNames(index);
    % Start forming the constraint string
    print_out = [constraintName{1},':'];
    for i=1:length(variables)
        if coefficients(1,i) == 1 && i ~= 1
            print_out = [print_out,' +'];
        elseif coefficients(1,i) == -1
            print_out = [print_out,' -'];
        elseif abs(coefficients(1,i)) ~= 1
            if coefficients(1,i) > 0
                print_out = [print_out,' +',num2str(coefficients(1,i))];
            elseif coefficients(1,i) < 0
                print_out = [print_out,' -',num2str(-coefficients(1,i))];
            end
        end
        print_out = [print_out,' ',variables{i,1}];
    end
    print_out = [print_out,' ',constraintType{index,1},' ',num2str(rhs(index,1))];
    if ScreenPrintflag
        fprintf([print_out,'\n'])
    end
    print_out_vec{j,1} = print_out;
end
