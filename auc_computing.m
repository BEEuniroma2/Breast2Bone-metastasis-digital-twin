function [AUC] = compute_auc(data,gt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
classes = unique(gt);
for f = 1 :size(data,2)

    [~,~,~,auc] = perfcurve(gt,data(:,f),classes(2));
    AUC(f) = auc; 
end

end