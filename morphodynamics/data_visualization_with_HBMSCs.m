% clear section
clear all; close all; clc;
addpath("violinplot_\");
%% load dataset
load("dataset_deep_features_with_HBMSCs.mat");
colors = [0.00,0.45,0.74;0.85,0.33,0.10];
%% compute auc of geometric features
id_sel = (GTE>=0);

auc_deep_mean = auc_computing(M_deep_mean(id_sel,:),GT_drug(id_sel));

figure()
plot(auc_deep_mean,"o")
axis([1,length(auc_deep_mean),0,1])
title("AUC deep features mean")




%% select features
rng(2)
th = 0.9;
id_sel_deep = (abs(auc_deep_mean-0.5)>quantile(abs(auc_deep_mean-0.5),th));
X = (M_deep_mean(id_sel,:)-mean(M_deep_mean(id_sel,:)))./std(M_deep_mean(id_sel,:));
Y = tsne(X(:,id_sel_deep));
[coeff,score,latent,tsquared,explained,mu] = pca(X(:,id_sel_deep));
%% plot
classes = unique(GT_drug(id_sel));
colours = colors;
marker = ["d","o"];
fontname = "Times New Roman";
titlefontsize = 20;
ylabelfontsize = 19;
axisfontsize = 12;
h(1) = figure();
for c = 1 : length(classes)
    id_class = GT_drug(id_sel)==classes(c);
    scatter(Y(id_class,2),Y(id_class,1),15,marker(c),"MarkerFaceColor",colours(c,:),"MarkerEdgeColor",colours(c,:))
    hold on
end
ax = gca;
ax.XAxis.FontSize = axisfontsize;
ax.YAxis.FontSize = axisfontsize;
ax.XAxis.FontName = fontname;
ax.YAxis.FontName = fontname;
xlabel("t-SNE variable 1","FontSize",ylabelfontsize,"FontName",fontname)
ylabel("t-SNE variable 2","FontSize",ylabelfontsize,"FontName",fontname)
title("t-SNE plot with HBMSCs","FontSize",titlefontsize,"FontName",fontname)

%%

m1 = mink(Y(:,1),10);
classes= unique(GT_drug);

id_control = (Y(:,1)==m1(1))&(GT_drug(id_sel)'==classes(1));
sum(id_control)
hold on 
plot(Y(id_control,2),Y(id_control,1),"*","MarkerSize",15,"Color",colours(1,:))
m2 = maxk(Y(:,1),10);

id_doxo = (Y(:,1)==m2(5))&(GT_drug(id_sel)'==classes(2));
plot(Y(id_doxo,2),Y(id_doxo,1),"*","MarkerSize",15,"Color",colours(2,:))
sum(id_doxo)

legend(["Control","Doxorubicin 50 nM"],"Location","northwest","FontSize",ylabelfontsize,"FontName",fontname)

gt_s = GTS(id_sel);
gt_c = GT_cell(id_sel);

id_subvideo_control = gt_s(id_control);
id_gt_s = find((gt_s==id_subvideo_control)&(GT_drug(id_sel)==classes(1)));



id_cell_control = find(id_control==1)-id_gt_s(1);

id_subvideo_doxo = gt_s(id_doxo);
id_gt_s_doxo = find((gt_s==id_subvideo_doxo)&(GT_drug(id_sel)==classes(2)));

id_cell_doxo = find(id_doxo==1)-id_gt_s_doxo(1);


base_dir = "D:\B2B\data analysis\3 motility and morphodynamics analysis\roi gallery deep\CTRL vs DOXO 50 nM with F37\";
dir_control = base_dir + num2str(id_subvideo_control)+"\"+num2str(id_cell_control,"%03d")+"\";

f_control = dir(dir_control+"*.png");

dir_doxo = base_dir + num2str(id_subvideo_doxo)+"\"+num2str(id_cell_doxo,"%03d")+"\";

f_doxo = dir(dir_doxo+"*.png");

image_control = imread(dir_control+f_control(1).name);

image_doxo = imread(dir_doxo+f_doxo(1).name);

f = figure();
imshow(image_control);
g = figure();
imshow(image_doxo);
savefig(f,"control_with_hbmcs_BC_cell.fig")

savefig(g,"doxo_with_hbmcs_BC_cell.fig")

imwrite(image_control,"control_with_hbmcs_BC_cell.png")
imwrite(image_doxo,"doxo_with_hbmcs_BC_cell.png")



