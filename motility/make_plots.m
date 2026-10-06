%% clear section
clear all; close all; clc;
addpath("violinplot_\")
%% load data
load("motility_dataset_all_classes_in_time.mat");
%% choose classes
% n = 1; % GelMA vs Laponite
% n = 2; % Control vs Doxorubicin without HBMSCs
n = 3; % Control vs Doxorubicin with HBMSCs

% time = 1; % 0h<t<12h
% time = 2; % 12h<t<24h
% time = 3; % 24h<t<36h
% time = 4; % 36h<t<48h
time = 5; % 48h<t<60h

if n == 1
    id_data = (GT_drug==0)&(GT_LAP == 1)&(GT_time==time)'; % look for the effect of hbmscs 
    y = GT_hbmscs(id_data);
    classes_names = ['CTRL', 'HBMSCS'];
    t = "GelMA+Laponite";
elseif n==2
    id_data = (GT_LAP == 1)&(GT_hbmscs==0)&(GT_time==time)'&((GT_drug==0)|(GT_drug==1)); % look for the effect of drug without hbmscs
    y = GT_drug(id_data);
    classes_names = {'Control', 'Doxorubcin 50 nM'};
    t = "without HBMSCs";
elseif n==3
    id_data = (GT_LAP == 1)&(GT_hbmscs==1)&(GT_time==time)'&((GT_drug==0)|(GT_drug==1)); % look for the effect of drug with hbmscs
    y = GT_drug(id_data);
    classes_names = {'Control', 'Doxorubcin 50 nM'};
    t = "with HBMSCs";

end


for i = 1 : length(y)
    GT{i} = classes_names{y(i)+1};
end

X = M(id_data,:);
y_time = GT_time(id_data);

res = 1.242*3; %um/px

f =1; % mean velocity
fontname = "Times New Roman";
titlefontsize = 20;
ylabelfontsize = 19;
axisfontsize = 12;
h(1) = figure();
violinplot(X(:,f)*res,GT');

ax = gca;
ax.XAxis.FontSize = ylabelfontsize;
ax.YAxis.FontSize = axisfontsize;
ax.XAxis.FontName = fontname;
ax.YAxis.FontName = fontname;
ylabel("cell velocity [µm/h]","FontSize",ylabelfontsize,"FontName",fontname)
if time == 1
    title_name = "0h < T < 12h";
elseif time ==2
    title_name = "12h < T < 24h";
elseif time ==3
    title_name = "24h < T < 36h";
elseif time ==4
    title_name = "36h < T < 48h";
elseif time == 5
    title_name = "Average cell velocity "+t;
end
title(title_name,"FontSize",titlefontsize,"FontName",fontname)


classes = unique(y);

[h2,p] = kstest2(X(y==classes(1),f),X(y==classes(2),f),"Tail","unequal");

