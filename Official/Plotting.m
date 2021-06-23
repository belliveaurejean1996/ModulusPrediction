%% Initialise programe
clear all %#ok<CLALL>
close all
clc

%% Load variable
load('D1_variable.mat')
b_D_1 = a;
a_D_1 = b;
Theta_D_1 = Phi;
Phi_D_1 = Theta;

load('D2_variable.mat')
b_D_2 = a;
a_D_2 = b;
Theta_D_2 = Phi;
Phi_D_2 = Theta;

%% Plot according Graphs
set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultAxesFontSize',30)

% figure(1)
% % histogram(b_D_2,'DisplayStyle', 'stairs')
% % hold on
% histogram(b_D_1,25)
% xlabel('Minor axis (pixels)')
% ylabel('Number of fibres')
% % legend
% % hold off
% 
% figure(2)
% % histogram(Phi_D_2)
% % hold on
% histogram(Phi_D_1,25)
% xlabel('\phi (pegrees)')
% ylabel('Number of fibres')
% % legend
% % hold off

figure(3)
histogram(Theta_D_1,30)
hold on
histogram(Theta_D_2,30)
xlabel('\theta (degrees)')
ylabel('Number of fibres')
legend('D1', 'D2')
hold off



