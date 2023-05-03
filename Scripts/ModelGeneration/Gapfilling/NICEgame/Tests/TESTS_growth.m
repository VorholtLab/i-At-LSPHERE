%% INPUT
%load Merged model with thermodynamics with defined media
load ('')
GFmodel = ;
%load the database used for gap-filling
load ('')
DBmodel = ;
% load 'WT' model with defined media
load('')
modelWT =;
% load the GF solutions
load('')
ActRxns = ;
changeCobraSolver('cplex_direct','LP');
%% find rxns from DB and block them
fromDB=(find(ismember(GFmodel.rxnFrom,{'DBmodel_1'})));

GFmodel.var_ub(ismember(GFmodel.varNames, strcat( 'F_',GFmodel.rxns(fromDB))))=0;
GFmodel.var_ub(ismember(GFmodel.varNames, strcat( 'R_',GFmodel.rxns(fromDB))))=0;
GFmodel.var_ub(ismember(GFmodel.varNames, strcat( 'NF_',GFmodel.rxns(fromDB))))=0;
GFmodel.var_lb(ismember(GFmodel.varNames, strcat( 'NF_',GFmodel.rxns(fromDB))))=0;


GFmodel.ub(fromDB)=0;
GFmodel.lb(fromDB)=0;
%% Get all alternatives in a list
 for i=1:size(ActRxns{1},1)
infocell{i,1} = concatenateList(ActRxns{1}{i},'|');
end
List=[infocell];
%% Test biomass and BBB yield FBA
[resultFBA]=TestFBA_growth(GFmodel,ActRxns{1});
%% Test biomass and BBB yield TFBA
[resultTFBA]=TestTFBA_growth(GFmodel,ActRxns{1});
%% Length
[resultLength]=TestLength_growth(ActRxns{1});         
%%  3 level EC for alternatives
% EC_WT is a mat file with the 3rd level EC numbers of the initial model
[resultEC]=TestEC_growth(GFmodel,DBmodel,ActRxns{1},EC_WT); 
%% Metabolites already in E.coli
[resultMets]=TestMets_growth(GFmodel,modelWT,ActRxns{1});
%% Table
resultFBA2=resultFBA(2:size(resultFBA,1),:);
resultTFBA2=resultTFBA(2:size(resultTFBA,1),:);
resultEC2=resultEC(2:size(resultEC,1),:);
resultMets2=resultMets(2:size(resultMets,1),:);
resultLength2=resultLength(2:size(resultLength,1),:);
table2{1,1}='Alternative Number';
table2{1,2}='Alternative';
table2{1,3}='Score FBA biomass';
table2{1,4}='Score TFBA biomass';
table2{1,5}='Score # of reactions';
table2{1,6}='Score EC number';
table2{1,7}='Score metabolites';
d=2;
for j=1:size(ActRxns{1},1)
table2{d,1}=j;
table2{d,2}=infocell{j,1};
table2{d,3}=resultFBA2{j,3};
table2{d,4}=resultTFBA2{j,3};
table2{d,5}=resultLength2{j,3};
table2{d,6}=resultEC2{j,5};
table2{d,7}=resultMets2{j,4};
d=d+1;
end

table2=cell2table(table2);