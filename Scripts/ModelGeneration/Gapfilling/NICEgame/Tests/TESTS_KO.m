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
rxnsToCheck={''};
ActRxns1={};
for i=1:length(rxnsToCheck)
    load(['gapFillingBiomass_' GFmodel.id '_' rxnsToCheck{i} '.mat'],'ActRxns')
ActRxns1=[ActRxns1; ActRxns] ;
end
ActRxns=ActRxns1;

changeCobraSolver('cplex_direct','LP');
%% find rxns from DB and block them

% DB KEGG or ATLAS
fromDB=(find(ismember(GFmodel.rxnFrom,{'DBmodel_1'})));

GFmodel.var_ub(ismember(GFmodel.varNames, strcat( 'F_',GFmodel.rxns(fromDB))))=0;
GFmodel.var_ub(ismember(GFmodel.varNames, strcat( 'R_',GFmodel.rxns(fromDB))))=0;

GFmodel.ub(fromDB)=0;
GFmodel.lb(fromDB)=0;
%% Check WT
% for KO gap-filing
changeCobraSolver('cplex_direct','LP');
%fba
solWTFBA=optimizeCbModel(GFmodel);
%thermo
solWT=solveTFAmodelCplex(GFmodel);
%% Get all alternatives in a list
names=rxnsToCheck';
for j=1:length(ActRxns)
    for i=1:size(ActRxns{j},1)
infocell{i,j} = concatenateList(ActRxns{j,1}{i},'|');
    end
end
List=[names; infocell];
%% Test biomass and BBB yield FBA
for i=1:length(rxnsToCheck)
[resultFBA]=TestFBA(GFmodel,rxnsToCheck(i),ActRxns{i},[],solWTFBA,[],TagBBB);
save(strcat('resultFBA_',rxnsToCheck{i},'.mat'),'resultFBA')
end
%% Test biomass and BBB yield TFBA
for i=1:length(rxnsToCheck)
[resultTFBA]=TestTFBA(GFmodel,rxnsToCheck(i),ActRxns{i},[],solWT,[],TagBBB);
save(strcat('resultTFBA_',rxnsToCheck{i},'.mat'),'resultTFBA')
end
%% Length
[resultLength]=TestLength(GFmodel,rxnsToCheck,ActRxns);                             
%%  3 level EC for alternatives
% EC_WT is a mat file with the 3rd level EC numbers of the initial model
[resultEC]=TestEC(GFmodel,modelWT,ATLAS,rxnsToCheck,ActRxns,EC_WT); 
%% Metabolites already in E.coli
[resultMets]=TestMets(GFmodel,modelWT,rxnsToCheck,ActRxns);
%% score
% get score
k=1;
for i=2:size(resultMets,1)
    if not(isempty(resultMets{i,6}))
        mets(k,1)=resultMets(i,6);
        k=k+1;
    end
end


k=1;
for i=2:size(resultEC,1)
    if not(isempty(resultEC{i,7}))
       EC(k,1)=resultEC(i,7);
        k=k+1;
    end
end
if TagBBB==0
    for i=2:size(resultFBA,1)
    resultFBA{i,102}=0;
    resultTFBA{i,102}=0;    
    end
end
     
jo=2;
for i=1:length(rxnsToCheck)
    k=1;
    for j=jo:jo+(size(ActRxns{i},1))-1
        score(k,i)=abs(resultFBA{j,6})+abs(resultFBA{j,102})+abs(resultTFBA{j,6})...
            +abs(resultTFBA{j,102})+abs(resultLength{j,4})+abs(EC{j-1,1})...
            +abs(mets{j-1,1});
        k=k+1;
    end
    jo=j+1;
end
for i=1:length(ActRxns)
l(i,1)=size(ActRxns{i},1);
end

rank1=zeros(max(l),length(ActRxns));
rank2=zeros(max(l),length(ActRxns));

for i=1:size(ActRxns,1)
    
rank11={};
rank22={};

[rank11,rank22]=sort(score(1:size(ActRxns{i},1),i));
rank1(1:size(ActRxns{i},1),i)=rank11;
rank2(1:size(ActRxns{i},1),i)=rank22;
end

d=1;
for i=1:size(rank1,2)
    rank(:,d)=rank2(:,i);
    rank(:,d+1)=rank1(:,i);
    d=d+2;
end
%% Table
resultFBA2=resultFBA(2:size(resultFBA,1),:);
resultTFBA2=resultTFBA(2:size(resultTFBA,1),:);
resultEC2=resultEC(2:size(resultEC,1),:);
resultMets2=resultMets(2:size(resultMets,1),:);
resultLength2=resultLength(2:size(resultLength,1),:);
d=1;
table2{1,1}='Reaction';
table2{1,2}='Alternative Number';
table2{1,3}='Alternative';
table2{1,4}='Score FBA biomass';
table2{1,5}='Score FBA BBBs';
table2{1,6}='Score TFBA biomass';
table2{1,7}='Score TFBA BBBs';
table2{1,8}='Score # of reactions';
table2{1,9}='Score EC number';
table2{1,10}='Score metabolites';
table2{1,11}='Total';
d=d+1;
jo=0;
for i=1:length(rxnsToCheck)
    for j=1:size(ActRxns{i},1)
table2{d,1}=rxnsToCheck{i};
table2{d,2}=rank2(j,i);
table2{d,3}=infocell{rank2(j,i),i};
table2{d,4}=resultFBA2{rank2(j,i)+jo,6};
table2{d,5}=resultFBA2{rank2(j,i)+jo,102};
table2{d,6}=resultTFBA2{rank2(j,i)+jo,6};
table2{d,7}=resultTFBA2{rank2(j,i)+jo,102};
table2{d,8}=resultLength2{rank2(j,i)+jo,4};
table2{d,9}=EC{rank2(j,i)+jo,1};
table2{d,10}=mets{rank2(j,i)+jo,1};
table2{d,11}=rank1(j,i);
d=d+1;
    end
    d=d+1;
    jo=jo+size(ActRxns{i},1);
end

table2=cell2table(table2);