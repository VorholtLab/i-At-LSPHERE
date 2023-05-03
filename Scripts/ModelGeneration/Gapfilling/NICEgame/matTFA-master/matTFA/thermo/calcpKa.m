function pKa_adj = calcpKa(name,index,ionicStr)
% calculate adjusted pKa value based on ionic strength
% name can be a specie name or pKa value itself

global pK;
global data;

if nargin == 2
   Kvalue = 1;
end

% if class(Kvalue) == 'double'
%    Kvalue = dbl2str(Kvalue);
% else
%    Kvalue = 1;
% end

% if name is a specie, then we check if it exists before going further

if existInStruct(pK,name)

    if class(name) == 'char'
        tmp = strcat('pK.',name,'{',dbl2str(index),'}');
        pKa = eval(tmp);
    end
    
    name = lower(name);
    charge = getspCharge(name,index);
    pH = 7;
    sigmanusq = 2*charge; %charge^2 + 1 - (charge-1)^2;
    lnkzero = log(10^-pKa);
    pKa_adj = -(lnkzero - (1.17582*ionicStr^0.5)*sigmanusq/(1+1.6*ionicStr^0.5))/log(10);

else

    pKa_adj = 'NA';
    fprintf(1,'specie %s not found in list of pKa values\n',name);
    return;

end
    

end