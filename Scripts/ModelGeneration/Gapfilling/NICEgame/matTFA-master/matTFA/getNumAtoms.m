function numAtoms = getNumAtoms(formula,type);

n = 1;
[compounds, tok] = regexp(formula, '([A-Z][a-z]*)(\d*)', 'match', 'tokens');

for j = 1:length(compounds) % go through each token.
    comp = tok{1,j}{1};
    
    if (length(tok{1,j}) == 1)
        q = 1;
    else
        q = str2num(tok{1,j}{2});
        
    end
    switch comp
        case 'C'
            Element{1,1} = 'C';
            Ematrix(n,1) = q;
        case 'N'
            Element{2,1} = 'N';
            Ematrix(n,2) = q;
        case 'O'
            Element{3,1} = 'O';
            Ematrix(n,3) = q;
        case 'H'
            Element{4,1} = 'H';
            Ematrix(n,4) = q;
        case 'P'
            Element{5,1} = 'P';
            Ematrix(n,5) = q;
        case 'S'
            Element{6,1} = 'S';
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Na'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Mg'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Cl'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'K'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Ca'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Mn'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Fe'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Ni'
            Ematrix(n,6) = Ematrix(n,6) + q;                
        case 'Co'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Cu'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Zn'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'As'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Se'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Ag'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Cd'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'W'
            Ematrix(n,6) = Ematrix(n,6) + q;
        case 'Hg'
            Ematrix(n,6) = Ematrix(n,6) + q;
%         otherwise
%             display('Warning');
%             display(formula)
%             display(comp);
    end
end

numAtoms = Ematrix(n,find(ismember(Element,type)));
