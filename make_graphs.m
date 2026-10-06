clear all; close all; clc;

%%
load("morphodynamics_dataset.mat");
colors = [0.00,0.45,0.74;0.85,0.33,0.10];
%%
id_sel = logical(ones(size(M_deep_no_proc_m,1),1));
M = M_deep_m;
auc = auc_computing(M(id_sel,:),gt_interaction(id_sel));
plot(auc,'o')
axis([1,length(auc),0,1])
%%

fontname = "Times New Roman";
titlefontsize = 22;
ylabelfontsize = 19;
axisfontsize = 12;


colours = colors;
marker = ["d","o","s","."];

id_sel = logical(ones(size(M,1),1));
th = 0.9;
id_sel_deep = (abs(auc-0.5)>quantile(abs(auc-0.5),th));



X = (M(id_sel,:)-mean(M(id_sel,:)))./std(M(id_sel,:));

[coeff,score,latent,tsquared,explained,mu] = pca(X(:,id_sel_deep));

%% 

y = gt_interaction;
classes = unique(y(id_sel));
colours = colors;

marker = ["d","o"];
figure()

v1 = 1;
v2 = 3;

for c = 1 : length(classes)
    id_class = y(id_sel)==classes(c);
    scatter(score(id_class,v1),score(id_class,v2),18,marker(c),"MarkerFaceColor",colours(c,:),"MarkerEdgeColor",colours(c,:))
    hold on
end
legend("Non interacting","Interacting","FontSize",ylabelfontsize-2,"FontName",fontname)
ax = gca;
ax.XAxis.FontSize = axisfontsize;
ax.YAxis.FontSize = axisfontsize;
ax.XAxis.FontName = fontname;
ax.YAxis.FontName = fontname;
xlabel("PC "+num2str(v1)+" ("+num2str(explained(v1),"%02.2f")+" %)","FontSize",ylabelfontsize,"FontName",fontname)
ylabel("PC "+num2str(v2)+ " ("+num2str(explained(v2),"%02.2f")+" %)","FontSize",ylabelfontsize,"FontName",fontname)
title("Morphodynamic features","FontSize",titlefontsize,"FontName",fontname)

%% 
figure()

for n =  1: length(gt_interaction)
    if gt_interaction(n) == 0
        GT{n} = 'Non interacting';
    elseif gt_interaction(n) == 1
        GT{n} = 'Interacting';
    end
end
X = mean_geom;

f = 5; % maj ax len


res = 1.242;
classes = unique(gt_interaction);
for i = 1 : length(classes)
    idx =  gt_interaction == classes(i);
    xdata = categorical(GT(idx))';
    ydata = X(idx, f) * res;
    boxchart(xdata, ydata,"BoxFaceColor",colors(i,:),'MarkerStyle','none');
    hold on
end
ax = gca;
ax.XAxis.FontSize = ylabelfontsize;
ax.YAxis.FontSize = axisfontsize;
ax.XAxis.FontName = fontname;
ax.YAxis.FontName = fontname;
title("Major axis length","FontSize",titlefontsize,"FontName",fontname)
ylabel("major axis length [μm]","FontSize",ylabelfontsize,"FontName",fontname)



[h, p] = kstest2(X(gt_interaction==classes(1),f),X(gt_interaction==classes(2),f),"Tail","unequal");
