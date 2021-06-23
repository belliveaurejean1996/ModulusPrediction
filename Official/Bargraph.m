% Data to be plotted as a bar graph
model_series = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
%Data to be plotted as the error bars
model_error = [1 2 3; 1 2 3; 1 2 3; 1 2 3];
% Creating axes and the bar graph
ax = axes;
h = bar(model_series,'BarWidth',1);
% Set color for each bar face
h(1).FaceColor = 'blue';
h(2).FaceColor = 'yellow';
% Properties of the bar graph as required
ax.YGrid = 'on';
ax.GridLineStyle = '-';
xticks(ax,[1 2 3 4]);
% Naming each of the bar groups
xticklabels(ax,{ 'Low', 'Middle', 'High', 'extra'});
% X and Y labels
xlabel ('Socio Economic Status');
ylabel ('Mean Writing Score');
% Creating a legend and placing it outside the bar plot
lg = legend('A','B','C','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end