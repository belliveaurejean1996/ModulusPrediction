%% Initialisation
clear all %#ok<CLALL>
clc
%#ok<*NOPTS>
%#ok<*SAGROW>

% data location
Start = 'D:\Université\Matrise\Article - Modulus Prediction\DataOutput\';
End = '_200x.xlsx';

% Name = ["R1","L1","R3","L3","RL5",...
%         "R1p","L1p","R3p","L3p","RT6",...
%         "BR1","BL1","L1L","L1R","L4L",...
%         "BR1p","BL1p","T2R_N2","T2L_N2","T3"];
Name = ["R1","L1","R3","RL5","RL5",...
        "R1p","L1p","R3p","L1p","L1p",...
        "BR1","BL1","L1L","L1R","BR1",...
        "T3","T2L_N2","T2R_N2","T2L_N2","T3"];


% Concatonate strings
files = strcat(Start,Name(1,:),End);

%% Material properties
% Ef = 220e9;     % Pa
Ef = 241e9;
E1 = 135e9;     % Pa
E2 = 10e9;      % Pa
G12 = 5.2e9;    % Pa
Volf = 0.59;    % Fiber volume fraction
Nu12 = 0.22;    % Estimation from composite book
Nu21 = (Nu12 * E2)/E1;
    
%% Facteur de correction
EA = 54.19e9;    % pa
EABD = 51.92e9;    % pa
E = 54e9;
FC1 = (EA/(Ef*Volf/3));
FC2 = (EABD/(Ef*Volf/3));
% FC = (E/(Ef*Volf/3))
FC = 1.21;

%% Modified laminate theorie
imax=length(Name);
percent=10;
for i = 1 : imax
    if i/imax>=percent/100
        disp([num2str(percent) ' %']);
        percent=percent+10;
    end
    
    % Read excel file
    Data = xlsread(files(1,i));

    % Calculate ellipse parameters
    [a,b,x,y,Phi,Theta,n] = EllipseParameters(Data);

    % A
    [Ex1_A] = ModifiedLaminate_A(E1,E2,G12,Volf,Nu12,Nu21,Phi,x,y,a,b);
    MethodeA(i) = Ex1_A/FC; 
%     MethodeA(i) = Ex1_A
    % ABD
    [Ex1_ABD] = ModifiedLaminate_ABD(E1,E2,G12,Volf,Nu12,Nu21,Phi,x,y,a,b);
    MethodeABD(i) = Ex1_ABD/FC;
%     MethodeABD(i) = Ex1_ABD
end

%% Model predictions
varnames = {'Direction','Samples','A','ABD'};
Direction = {'RD1','RD1','RD1','RD1','RD1',...
             'RD2','RD2','RD2','RD2','RD2',...
             'BD1','BD1','BD1','BD1','BD1',...
             'BD2','BD2','BD2','BD2','BD2'};
Results = table(Direction.',Name.',MethodeA.', MethodeABD.','VariableNames',varnames);

% Means
MAL = mean(MethodeA(1:5));
MAS = mean(MethodeA(6:10));
MABL = mean(MethodeA(11:15));
MABS = mean(MethodeA(16:20));

MABDL = mean(MethodeABD(1:5));
MABDS = mean(MethodeABD(6:10));
MABDBL = mean(MethodeABD(11:15));
MABDBS = mean(MethodeABD(16:20));

% STD
SAL = std(MethodeA(1:5));
SAS = std(MethodeA(6:10));
SABL = std(MethodeA(11:15));
SABS = std(MethodeA(16:20));

SABDL = std(MethodeABD(1:5));
SABDS = std(MethodeABD(6:10));
SABDBL = std(MethodeABD(11:15));
SABDBS = std(MethodeABD(16:20));


%% load experimental data
Exp = xlsread('D:\Université\Matrise\Article - Modulus Prediction\DataOutput\Experimental.xlsx');

ExpLong = Exp(:,1);
ExpShort = Exp(:,2);
ExpBaseLong = Exp(:,3);
ExpBaseShort = Exp(:,4);

MEL = mean(ExpLong);
MES = mean(ExpShort);
MEBL = mean(ExpBaseLong);
MEBS = mean(ExpBaseShort);

SEL = std(ExpLong);
SES = std(ExpShort);
SEBL = std(ExpBaseLong);
SEBS = std(ExpBaseShort);


%% figures
figure('DefaultAxesFontSize',16)
% bar chart data
% BarMeans = [MEBL,MABL,MABDBL; MEBS,MABS,MABDBS; MEL,MAL,MABDL; MES,MAS,MABDS];
% BarSTD = [SEBL,SABL,SABDBL; SEBS,SABS,SABDBS; SEL,SAL,SABDL; SES,SAS,SABDS];

BarMeans = [MEBL,MABDBL; MEBS,MABDBS; MEL,MABDL; MES,MABDS];
BarSTD = [SEBL,SABDBL; SEBS,SABDBS; SEL,SABDL; SES,SABDS];

% set(0,'defaultAxesFontSize',24)
% Creating axes and the bar graph
ax = axes;
h = bar(BarMeans,'BarWidth',1);

% Set color for each bar face
h(1).FaceColor = [0,0,0];
h(2).FaceColor = [0.75 0.75 0.75];
% h(3).FaceColor = [1,1,1];

% Properties of the bar graph as required
xticks(ax,[1 2 3 4]);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',20)
set(gca,'XTickLabelMode','auto')
% Naming each of the bar groups
xticklabels(ax,{'Baseline [D1]','Baseline [D2]','Recycled [D1]','Recycled [D2]'});


% X and Y labels
ylabel('Young''s Modulus [GPa]','FontSize',24);

% Creating a legend and placing it outside the bar plot
lg = legend('Mechanical testing','Model Prediction','AutoUpdate','off','FontSize',20);
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;

% Finding the number of groups and the number of bars in each group
ngroups = size(BarMeans, 1);
nbars = size(BarMeans, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, BarMeans(:,i), BarSTD(:,i), 'Color', [0.35 0.35 0.35], 'linestyle', 'none');
end

% reference lines
yline(52,'--','QI reference','LabelHorizontalAlignment','left','FontSize',20,'LineWidth',2);
% yline(43,'--','MC 1200','LabelHorizontalAlignment','left','FontSize',18);
hold off

%% statistical analysis

RD1_A = MethodeA(1:5).';
RD2_A = MethodeA(6:10).';
BD1_A = MethodeA(11:15).';
BD2_A = MethodeA(16:20).';

RD1_ABD = MethodeABD(1:5).';
RD2_ABD = MethodeABD(6:10).';
BD1_ABD = MethodeABD(11:15).';
BD2_ABD = MethodeABD(16:20).';

% anova analysis
exp = repmat("expirimental",[length(ExpLong), 1]);
A = repmat("A",[length(RD1_A), 1]);
ABD = repmat("ABD",[length(RD1_ABD), 1]);

% Baseline D1
Groupe = [exp; A];
BD1 = [ExpBaseLong; BD1_ABD];
P_BD1 = anova1(BD1, Groupe,'off');

% Baseline D2
Groupe = [exp; A];
BD2 = [ExpBaseShort; BD2_ABD];
P_BD2 = anova1(BD2, Groupe,'off');

% Recycled D1
Groupe = [exp; A];
RD1 = [ExpLong; RD1_ABD];
P_RD1 = anova1(RD1, Groupe,'off');

% Recycled D1
Groupe = [exp; A];
RD2 = [ExpShort; RD2_ABD];
P_RD2 = anova1(RD2, Groupe,'off');

% Results anova
P_value = [P_BD1; P_BD2; P_RD1; P_RD2]; 
Samples = ["BD1";"BD2";"RD1";"RD2"];

Alpha = 0.05;
for i = 1:length(P_value)
    if P_value(i) > Alpha
        Difference(i,1) = "no";
    else
        Difference(i,1) = "Yes";
    end
end

Stats = table(Samples, P_value, Difference)








