function [AllTransports, TransportNoCouples, CoupledTransports, ImportantTransports, directions, TransportGroups] = identifyTransportRxns(model,biomassRxnNames,ATPsynth_RxnNames)

% this function identifies all the transport reactions in a model and sorts
% them according to compartments
% Milenko Tokic 02.11.2015

% INPUT:
%   - cobra frendly model, metabolites should be stored in a field model.mets
%     and field model.metCompSymbol is also required

% OUTPUT
%   - AllTransports: list of a transport reactions, formulas and compartments between
%     metabolite is being exchanged

%   - TransportNoCouples: table with list of transport reaction without
%     carbon couples (e.g. agm_c + arg-L_p <=> arg-L_c + agm_p )

%   - List of carbon coupled transports (e.g. agm_c + arg-L_p <=> arg-L_c + agm_p)

%   - ImportantTransports: list of transport that we have to assign
%     directionality

%   - directions: vector of 1 and -1, based on this we will assign
%     directionality

%   - TransportGroups: reactions grouped into bims according compartment
%     info

% remove biomass reaction, this is not a transport and we dont want to keep
% it becouse it will create a mess

for i = 1:length(biomassRxnNames)
    indexBiomass = find(~cellfun(@isempty, regexp(model.rxns, biomassRxnNames{i})));
    if isempty(indexBiomass)==0
        model.rxns(indexBiomass) = [];
        model.S(:,indexBiomass) = [];
    end
end
%%
% first we remove compartment information

metsNoComp = cellfun(@(x)(x(1:end-2)), model.mets, 'uni', false);

% find the list of metabolites in each reaction
% first we want to identify nonzero elements in a reaction

NonZero = {};
for i = 1:length(model.rxns)
    INDEX = find(model.S(:,(i)));
    NonZero{i} = INDEX;
end

% put metabolite names instead of indeces

MetsInRxns = {};
for i = 1:length(NonZero)
    indentifier = NonZero{i};
    MetsInRxns{i} = metsNoComp(indentifier);
end


% find unique values in a MetsInRxns
Unique = {};
for i = 1:length(MetsInRxns)
    Unique{i} = unique(MetsInRxns{:,i},'stable');
end
%%
% Now, we want also to group reaction based on compartments. For example if the
% metabolite is being transported from extracellular space to periplasm, we
% will store it in EP group etc.

% First, we want to extract compartment data

metsComp = cellfun(@(x)(x(end)),model.mets, 'uni', false);

% now we want to put compartment info in a reaction
MetsCompInRxns = {};
for i = 1:length(NonZero)
    indentifier = NonZero{i};
    MetsCompInRxns{i} = metsComp(indentifier);
end

% find unique values in a MetsCompInRxns
UniqueComp={};
for i = 1:length(MetsInRxns)
    UniqueComp{i} = unique( MetsCompInRxns{:,i});
end

% identify transport reactions by comparing Uniqe and MetsInRxns.
% if there is a difference, that means that at least one metabolite
% participates from both side of the equations, meaning, it is being
% transported

Transport_reactions = [];
for i = 1:length(MetsInRxns)
    if size(MetsInRxns{i},1)==size(Unique{i},1)
        % Reaction is not a transport, we store 0
        Transport_reactions(i) = 0;
    else
        % Transport reaction! we store 1,
        Transport_reactions(i) = 1;
    end
    
end

% identify transport reactions and their copartment info
Transport_reactions = Transport_reactions';
IndexTransport = find(Transport_reactions);
AllTransports = model.rxns(IndexTransport);
compInfo = UniqueComp(find(Transport_reactions));
compInfo = compInfo';
%%
% Exclute ATP synthase from the transports.

for i = 1:length(ATPsynth_RxnNames)
    indexATPS = find(~cellfun(@isempty, regexp(AllTransports, ATPsynth_RxnNames{i})));
    AllTransports(indexATPS,:) = [];
    IndexTransport(indexATPS,:) = [];
    compInfo(indexATPS,:) = [];
end


% merge compartment data

for i=1:length(compInfo)
    compInfo{i,1} = compInfo{i,1}';
    compInfo{i,1} = strjoinTGL(compInfo{i,1});
    compInfo{i,1} = strrep(compInfo{i,1},' ','');
end
%%
% some ractions can be identified as transport throuh 3 compartmets like:
% K2L4Aabctex	atp_c + h2o_c + kdo2lipid4_p <=> adp_c + h_c + pi_c + kdo2lipid4_e
% this reaction will apear in the cpe group because of the h, so the info
% for h will be removed

for i=1:length(compInfo)
    MultiComp(i) = length(cell2str(compInfo(i)));
    IndexToChangeA = find(MultiComp>2);
end

% this part will be executed only if there is such a case (three
% compartments)

if isempty(IndexToChangeA)==0
    TransportToChange = AllTransports(IndexToChangeA);
    
    for i = 1:length(TransportToChange)
        [~,IndexToChange(i)] = ismember(TransportToChange(i),model.rxns);
    end
    
    % identify unique metabolite for these reactions and store its
    % compartment
    
    for i = 1:length(IndexToChange)
        
        IDStriple{i} = find(model.S(:,IndexToChange(i)) < 0);
        IDPtriple{i} = find(model.S(:,IndexToChange(i)) > 0);
        
    end
    
    for i = 1:length(IDStriple)
        SubTriple{:,i}  = metsNoComp(IDStriple{:,i});
        ProdTriple{:,i} = metsNoComp(IDPtriple{:,i});
    end
    
    for i=1:length(SubTriple)
        TranspotedMetabolitetriple{i} = intersect(SubTriple{:,i},ProdTriple{:,i});
    end
    
    % delete hydrogen from these and now find original compartment
    
    for i = 1:length(TranspotedMetabolitetriple)
        [m(i) d(i)] = size(TranspotedMetabolitetriple{:,i});
        if m(i)>1
            TranspotedMetabolitetriple{:,i} = setdiff(TranspotedMetabolitetriple{:,i},'h');
        else
            TranspotedMetabolitetriple{:,i} = TranspotedMetabolitetriple{:,i};
        end
    end
    
    % restore the original compartment
    
    for i=1:length(IDStriple)
        SubTripleReal{:,i} = model.mets(IDStriple{:,i});
        ProdTripleReal{:,i} = model.mets(IDPtriple{:,i});
        participatnts{i} = union(SubTripleReal{:,i},ProdTripleReal{:,i});
    end
    
    for i = 1:length(TranspotedMetabolitetriple)
        indexMetsToChange{i} = find(~cellfun(@isempty, regexp(participatnts{i}, TranspotedMetabolitetriple{i})));
    end
    
    for i = 1:length(participatnts)
        N = participatnts{i};
        for j = indexMetsToChange{i}
            UniqueParticipants{i} = N(j);
        end
    end
    
    for i = 1:length(UniqueParticipants)
        CorComInfo{i} = cellfun(@(x)(x(end)),UniqueParticipants{i}, 'uni', false);
    end
    
    % put everything in one cell
    
    for i = 1:size((CorComInfo),2)
        CorComInfo1 = (CorComInfo{:,i});
        Mergecells{i,1} = strjoinTGL(CorComInfo1');
        compInfoToChange{i,1} = strrep(Mergecells{i,1},' ','');
    end
    
    % find Indexes in compInfo to be changed
    
    for i = 1:length(compInfoToChange)
        compInfo{IndexToChangeA(i)} = compInfoToChange{i};
    end
end

%%

% now, we want to indentify unique metabolite for all the other reactions that is being transported
% identify substrates and products in transport reaction

% find indexes of transport reactions


for i = 1:length(IndexTransport)
    
    IDS{i} = find(model.S(:,IndexTransport(i)) < 0);
    IDP{i} = find(model.S(:,IndexTransport(i)) > 0);
    
end

% put metabolite names instead of indeces

for i = 1:length(IDS)
    Sub{:,i}  = metsNoComp(IDS{:,i});
    Prod{:,i} = metsNoComp(IDP{:,i});
end

for i=1:length(Sub)
    TranspotedMetabolite{i} = intersect(Sub{:,i},Prod{:,i});
end

% if there is more than two elements in the cell that means that at least
% one is h and it should be removed, also empty cells

for i = 1:length(TranspotedMetabolite)
    [m(i) d(i)] = size(TranspotedMetabolite{:,i});
    if m(i)>1
        TranspotedMetabolite1{:,i} = setdiff(TranspotedMetabolite{:,i},'h');
    else
        TranspotedMetabolite1{:,i} = TranspotedMetabolite{:,i};
    end
end

% second itteration, removing na1

for i = 1:length(TranspotedMetabolite1)
    [m(i) d(i)] = size(TranspotedMetabolite1{:,i});
    
    if m(i)>1
        TranspotedMetabolite2{:,i} = setdiff(TranspotedMetabolite1{:,i},'na1');
    else
        TranspotedMetabolite2{:,i} = TranspotedMetabolite1{:,i};
    end
end

% third itteration, removing pi

for i = 1:length(TranspotedMetabolite2)
    [m(i) d(i)] = size(TranspotedMetabolite2{:,i});
    
    if m(i)>1
        TranspotedMetaboliteUnique{:,i} = setdiff(TranspotedMetabolite2{:,i},'pi');
    else
        TranspotedMetaboliteUnique{:,i} = TranspotedMetabolite2{:,i};
    end
end


TranspotedMetaboliteUnique = TranspotedMetaboliteUnique';


TranspotedMetaboliteUnique(indexATPS) = [];
% delete empty cells
emptyCells = cellfun(@isempty,TranspotedMetaboliteUnique);
emptyCells = find(emptyCells);
TranspotedMetaboliteUnique(emptyCells) = [];



AllTransports(emptyCells) = [];
compInfo(emptyCells) = [];

% print Rxn Formula for visual inspection

Transports_Formulas = printRxnFormula(model,AllTransports);

% storing everything in a onr matrix

AllTransports = [AllTransports Transports_Formulas TranspotedMetaboliteUnique compInfo];

%% remove carbon couples (e.g. agm_c + arg-L_p <=> arg-L_c + agm_p )

Metab1 = AllTransports(:,3);

for i = 1:length(Metab1)
    A(i) = size((Metab1{i}),1);
end
IdentifyAntiports = find(A>1)';
TransportNoCouples = AllTransports;

TransportNoCouples(IdentifyAntiports,:) = [];

Metab2 = TransportNoCouples(:,3);
for i = 1:length(Metab2)
    MetabMat{i} = cell2mat(Metab2{i});
end
TransportNoCouples(:,3) = MetabMat';
CoupledTransports = AllTransports(IdentifyAntiports,:);


%% Grouping by compartments

if isfield(model,'metCompSymbol')==1
    transportGroups = unique(model.metCompSymbol');
else
    transportGroups = unique(metsComp);
end

transportGroups = combntns(transportGroups,2);
for i = 1:length(transportGroups)
    transportGroupsJoined{i,1} = strcat(transportGroups{i,1}, transportGroups{i,2});
end

for i = 1:length(AllTransports)
    for j = 1:length(transportGroupsJoined)
        Index = strncmp(AllTransports{i, 4}, transportGroupsJoined{j}, 3);
        if Index==1
            BIM{i, j} = AllTransports{i, :};
        else
            BIM{i, j} = [];
        end
    end
end

% remove empty cells from the table
% first, separte them into different columns

Names = transportGroupsJoined;

for i = 1 : size(BIM(1,:),2)
    eval(['A',num2str(i),' = BIM(:, i);']);
    eval(['BIMS.A',num2str(i),' = A',num2str(i),'(~cellfun(@isempty,A',num2str(i),'));']);
end

% put names of compartments instead the generic once

for i = 1:length(Names)
    eval(['TransportGroups.' ,Names{i} ,' =BIMS.A',num2str(i),';']);
end

%% Now we want to indentify transports that we have to assign directionality

% if a metabolite in the same bim appears more than ones, that means that
% the directionality of this transports has to be assigned

Futile = {};
CompBim = unique(TransportNoCouples(:, 4));
for i = 1: length(CompBim)
    A = ismember(TransportNoCouples(:,4), CompBim{i});
    B = find(A==1);
    C = TransportNoCouples([B],[1,3,4]);
    ListMetabolites = C(:,2);
    ListReactions = C(:,1);
    ListCompartments = C(:,3);
    
    ListMetabolitesUnique = unique(ListMetabolites);
    if size(ListMetabolitesUnique,1) == size(ListMetabolites,1)
        Futile{i} = {};
    else
        [position] = findDoubles(ListMetabolites);
        RepetedMetabolites = ListMetabolites(position);
        RepetedMetabolitesRxns = ListReactions(position);
        RepetedMetabolitescompartments = ListCompartments(position);
        Futile{i} = [RepetedMetabolitesRxns,RepetedMetabolites,RepetedMetabolitescompartments];
    end
end
ImportantTransportsTemp = {};

for i=1:length(Futile)
    eval(['importantTransports',num2str(i),' = Futile{:, i};']);
    CurTime = eval(horzcat('importantTransports', num2str(i)));
    ImportantTransportsTemp = [ImportantTransportsTemp; CurTime];
end


% remove h from important Transports
[~,hpos] = ismember(ImportantTransportsTemp(:,2), 'h');
hpos = find(hpos==1);
ImportantTransportsTemp(hpos, :) = [];
ImportantTransports = sortrows(ImportantTransportsTemp);


% Check the soichiometry, find direction of an important transport couple
% and put the same number if they operate in the same direction.
% If the numbers are -1 and 1 that means that they do not operate in the same direction.
% in this case 1 means that the direction is E -> P -> C, and -1 opposite.

for i = 1:size(ImportantTransports(:,1),1)
    [~, IndexToChange(i)] = ismember(ImportantTransports(i,1),model.rxns);
end
IndexToChange = IndexToChange';
directions = [];
A = full(model.S);
for i= 1:length(IndexToChange)
    participants{i} = find(A(:, IndexToChange(i)));
    Indexes = A(participants{i}, IndexToChange(i));
    Mets = model.mets(find(A(:, IndexToChange(i))));
    comp = [Mets, num2cell(Indexes)];
    sorted = sortrows(comp, 2);
    idx = strfind(sorted(:, 1), cell2mat(ImportantTransports(i,2)));
    IndNonEmpty = find(~cellfun(@isempty, idx));
    TransMets = sorted((IndNonEmpty),1);
    metsNoComp = cellfun(@(x)(x(end)), TransMets, 'uni', false);
    
    metsNoComp=metsNoComp';
    metsNoComp = strjoinTGL(metsNoComp);
    metsNoComp = strrep(metsNoComp,' ','');
    
    if  isequal(metsNoComp,ImportantTransports{i,3})==0
        directions(i) = 1;
    else
        directions(i) = -1;
    end
end

ImportantTransports = [ImportantTransports,num2cell(directions')];

end
