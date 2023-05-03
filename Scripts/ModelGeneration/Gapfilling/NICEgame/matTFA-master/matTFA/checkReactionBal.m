function [model,result] = checkReactionBal(model,rxns,addProtons)
% checks if reaction is balanced
% basic checks

sum_charge = 0;
sum_charge_noH = 0;

rxn_index = find(ismember(model.rxns,rxns));

if length(rxn_index) > 1
    fprintf('duplicate reaction name found:\n');
    disp(rxns)
    keyboard
end

mets = model.mets(find(model.S(:,rxn_index)));
metFormulas = model.metFormulas(find(model.S(:,rxn_index)));
metCharges = model.metCharge(find(model.S(:,rxn_index)));
stoich = model.S(find(model.S(:,rxn_index)),rxn_index);
Ematrix = zeros(length(mets), 6);

if length(mets) == 1    
    result = 'drain flux';
    
elseif ~isempty(find(ismember(metFormulas,'NA')))
    result = 'missing structures';
    return
    
elseif length(mets) == length(stoich)

    for n=1:length(mets)
        met_index = find(ismember(model.mets,mets(n)));
        if length(met_index > 1)
            met_index = met_index(1);
        end
        
        if isempty(met_index)
        	error('metabolite not found');
        end
        
        met_formula = model.metFormulas(met_index);
        met_charge = model.metCharge(met_index);
        
        if isempty(met_formula) || isempty(met_charge)
        	error('metabolite charge or formula not found');
        end
        
        if ~strcmp(met_formula,'H')
            sum_charge_noH = sum_charge_noH + stoich(n)*met_charge;
        end
        
        sum_charge = sum_charge + stoich(n)*met_charge;
        
        [compounds, tok] = regexp(met_formula, '([A-Z][a-z]*)(\d*)', 'match', 'tokens');
        tok = tok{1,1};
        
        for j = 1:length(tok) % go through each token.
            t = tok{1,j};
            comp = t{1,1};
            q = str2num(t{1,2});
            if (isempty(q))
                q = 1;
            end
            
            switch comp
            case 'C'
                Ematrix(n,1) = q;
            case 'N'
                Ematrix(n,2) = q;
            case 'O'
                Ematrix(n,3) = q;
            case 'H'
                Ematrix(n,4) = q;
            case 'P'
                Ematrix(n,5) = q;
            case 'Na'
                Ematrix(n,6) = q;
            case 'Mg'
                Ematrix(n,7) = q;
            case 'S'
                Ematrix(n,8) = q;
            case 'Cl'
                Ematrix(n,9) = q;
            case 'K'
                Ematrix(n,10) = q;
            case 'Ca'
                Ematrix(n,11) = q;
            case 'Mn'
                Ematrix(n,12) = q;
            case 'Fe'
                Ematrix(n,13) = q;
            case 'Ni'
                Ematrix(n,14) = q;                
            case 'Co'
                Ematrix(n,15) = q;
            case 'Cu'
                Ematrix(n,16) = q;
            case 'Zn'
                Ematrix(n,17) = q;
            case 'As'
                Ematrix(n,18) = q;
            case 'Se'
                Ematrix(n,19) = q;
            case 'Ag'
                Ematrix(n,20) = q;
            case 'Cd'
                Ematrix(n,21) = q;
            case 'W'
                Ematrix(n,22) = q;
            case 'Hg'
                Ematrix(n,23) = q;
            case 'R'
                Ematrix(n,24) = q;
            case 'Mo'
                Ematrix(n,25) = q;
            case 'X'
                Ematrix(n,26) = q;
            otherwise
%                 display('Warning');
%                 display(met_formula)
%                 display(compounds);
            end
        end
    end
    
    bal = stoich'*Ematrix;

    if isempty(find(bal)) && (sum_charge == 0)
        result = 'balanced';
    elseif (length(find(bal)) == 1)
        if (find(bal) == 4) && (bal(4) == sum_charge)
            rxn_index = find(ismember(model.rxns,rxns));
            rxn_comp = model.rxnComp(rxn_index);
            
            % search using metSEEDID and metCompSymbol is more accurate
            % otherwise we try to search by formula
            if isfield(model,'metSEEDID')
                proton_met_index = find(ismember(model.metSEEDID,'cpd00067'));
            elseif isfield(model,'metFormula')
                proton_met_index = find(ismember(model.metFormula,'H'));
                if length(proton_met_index) ~= 1
                    fprintf('Cannot find proton metabolite index\n');
                end
            else
                fprintf('Cannot identify proton\n');
            end
            
            if isfield(model,'metCompSymbol')
                proton_index = proton_met_index(find(ismember(model.metCompSymbol(proton_met_index),rxn_comp)));
            else
                fprintf('Cannot find compartment info\n');
            end

            if isempty(proton_index)
                result = 'missing atoms';
                return
            elseif (addProtons)
                if (sum_charge_noH > 0)
                    model.S(proton_index,rxn_index) = -sum_charge_noH;
                    result = [int2str(sum_charge_noH) ' protons added to reactants'];
                else
                    model.S(proton_index,rxn_index) = -sum_charge_noH;
                    result = [int2str(-sum_charge_noH) ' protons added to products'];
                end
            else
                result = 'protons not added';
            end
        else
            result = 'missing atoms';
        end
        
    else
        result = 'missing atoms';
    end
end
    
    
    
        