function [FMeasure, Recall, Precision, RI, FMI] = compute_matching_matrix(U0, U)
% ����U0��U������FMeasure, Recall, Precision, RI, FMI��ֵ
% U0Ϊ�ο����־���
% U0��U��Ϊcxn����
% ref:
% 1. Michele Ceccarelli, Antonio Maratea. Improving fuzzy clustering of
% biological data by metric learning with side information. 
% 2. Jesse Davis, Mark Goadrich. The Relationship Between Precision-Recall
% and ROC Curves
% last modified: 2009.7.6 18:54:00

[cluster_n0, data_n0] = size(U0);
[cluster_n, data_n] = size(U);

if (data_n0~=data_n)
    fprintf('FMeasure: �������С��ͬ\n');
    FMeasure = 0;
    return;
end

if cluster_n < cluster_n0
    zeroM = zeros(cluster_n0-cluster_n, data_n);
    U = [U; zeroM];
else
    if cluster_n > cluster_n0
        zeroM = zeros(cluster_n-cluster_n0, data_n);
        U0 = [U0; zeroM];
    end
end

mask = triu(ones(data_n0, data_n0))-eye(data_n0);
% ����Mr
M0 = (U0'*U0).*mask;
% ����M
M = (U'*U).*mask;
% ����TP,FP,FN,TN


TP = sum(sum(M0 & M));
FN = sum(sum(M0 > M));
FP = sum(sum(M0 < M));
TN = sum(sum(1-M0 & 1-M))-(data_n0+1)*data_n0/2;

Recall = TP/(TP+FN);
Precision = TP/(TP+FP);
FMeasure = 2*Recall*Precision/(Recall+Precision);
RI = (TP+TN)/(TP+TN+FP+FN);
FMI = sqrt((TP/(TP+FP))*(TP/(TP+FN)));
