clear all; close all;
clc;
load("dataset_motility_interaction.mat");
colors = [0.00,0.45,0.74;0.85,0.33,0.10];

classes = unique(GT_drug);
for n_class = 1 : length(classes)
    id_class = (GT_drug == classes(n_class));
    classes_ = unique(GT_interaction(id_class));
    for n = 1 : length(classes_)

        h = find((GT_interaction==classes_(n))&(id_class)==1);
        for j = 1 : length(h)
            if n_class == 1 && n == 1
                GT{h(j)} = 'Control non interacting';
                gt(h(j)) = 1;
            elseif n_class == 1 && n ==2
                GT{h(j)} = 'Control interacting';
                gt(h(j)) = 3;
            elseif n_class == 2 && n == 1

                GT{h(j)} = 'Doxorubicin non interacting';
                gt(h(j)) = 2;
            elseif n_class == 2 && n == 2
                GT{h(j)} = 'Doxorubicin interacting';
                gt(h(j)) = 4;
            end
        end
    end
end
m = median(M(gt==1,1));


fontname = "Times New Roman";
titlefontsize = 22;
ylabelfontsize = 19;
axisfontsize = 12;
res = 1.242*3;%um/px
%%




%%
for n = 1 : length(GT_interaction)
    if GT_interaction(n) == 0
        GT2{n} = 'Non interacting';
    elseif GT_interaction(n) == 1
        GT2{n} = 'Interacting';
    end
end

res = 1.242*3;%um/px
m = median(M(gt==1,1))*res;
% id_sel = GT_drug==0;
f = 1;
% classes = unique(gt);
classes = unique(GT_interaction);

for i = 1 : length(classes)
    if i == 1
        c = [0 0.3470 0.6410];
    elseif i==3
        c = [0 0.4470 0.9];
    elseif i == 4
        c = [0.7500 0.2250 0.0980];
    elseif i ==2
        c= [0.980 0.3250 0.0980];
    end
    boxchart(categorical(GT2(GT_interaction==classes(i))), M(GT_interaction==classes(i),f)*res,"BoxFaceColor",colors(i,:),'MarkerStyle','none');
    hold on
end
ax = gca;
ax.XAxis.FontSize = ylabelfontsize;
ax.YAxis.FontSize = axisfontsize;
ax.XAxis.FontName = fontname;
ax.YAxis.FontName = fontname;
ylabel("Avg. approaching velocity [μm/h]","FontSize",ylabelfontsize-1,"FontName",fontname)
title("Average approaching velocity","FontSize",titlefontsize,"FontName",fontname)

[h,p] = kstest2(M(GT_interaction==classes(1),f),M(GT_interaction==classes(2),f),"Tail", "larger");

%%
figure()
f = 9;
classes = unique(GT_interaction);
res = 1/1.242;
for i = 1 : length(classes)
    if i == 1
        c = [0 0.3470 0.6410];
    elseif i==3
        c = [0 0.4470 0.9];
    elseif i == 4
        c = [0.7500 0.2250 0.0980];
    elseif i ==2
        c= [0.980 0.3250 0.0980];
    end
    boxchart(categorical(GT2(GT_interaction==classes(i))), M(GT_interaction==classes(i),f),"BoxFaceColor",colors(i,:),'MarkerStyle','none');
    hold on
end
ylabel("Curvature [1/μm]")
title("Curvature")
ax = gca;
ax.XAxis.FontSize = ylabelfontsize;
ax.YAxis.FontSize = axisfontsize;
ax.XAxis.FontName = fontname;
ax.YAxis.FontName = fontname;
ylabel("Curvature [1/μm]","FontSize",ylabelfontsize,"FontName",fontname)
title("Curvature","FontSize",titlefontsize,"FontName",fontname)


ylim([0,4])

[h,p] = kstest2(M(GT_interaction==classes(1),f),M(GT_interaction==classes(2),f),"Tail","smaller");


