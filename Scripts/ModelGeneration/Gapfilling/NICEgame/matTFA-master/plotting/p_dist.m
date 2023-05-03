function p_dist(x)

BINCOUNT = 50;
LEGEND_LOC = 'southwest';

% Fits for side histograms
[hist1, edges1] = histcounts(x,BINCOUNT);
pd_1_n = fitdist(x,'normal');
pd_1_k = fitdist(x,'kernel');
centers1 = (edges1(1:end-1)+edges1(2:end))/2;

% Create figure
figure1 = figure;
figure1.Position = [50,800,1500,500];
% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.336117021276596 0.11 0.362765957446809 0.815]);
hold(axes1,'on');


ls = linspace(min(x), max(x));

% Create subplot
subplot1 = subplot(1,1,1,'Parent',figure1);
hold(subplot1,'on');

% Create bar
bar(centers1,hist1,'DisplayName','residual count','Parent',subplot1,...
    'FaceColor',[0 0.447058826684952 0.74117648601532],...
    'EdgeColor',[0.0784313753247261 0.168627455830574 0.549019634723663],...
    'BarWidth',1);
% Define plot data
dx1 = edges1(2)-edges1(1);

dfit_1_n = pdf(pd_1_n,ls)*length(x)*dx1;
dfit_1_k = pdf(pd_1_k,ls)*length(x)*dx1;

% Create kernel plot
plot(ls,dfit_1_k,'Parent',subplot1,'DisplayName','fitted residual distribution',...
    'LineWidth',2,...
    'Color',[0 0 0]);
% % Create normal plot
% plot(ls,dfit_1_n,'Parent',subplot1,'DisplayName','expected normal distribution',...
%     'LineWidth',2,...
%     'LineStyle','--',...
%     'Color',[1 0 0]);

