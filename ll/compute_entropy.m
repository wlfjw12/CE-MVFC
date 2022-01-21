function Entropy = compute_entropy(U0, U)
% ����U0��U������Entropy
% U0Ϊ�ο����־���
% �ο����ף�Liping Jing, Michael K. Ng, Joshua Zhexue Huang
% An Entropy Weighting k-Means Algorithm for Subspace Clustering of High-Dimensional Sparse Data. 
% IEEE Trans. Knowl. Data Eng. (TKDE) 19(8):1026-1041 (2007)
% Updated: 2009.07.06 19:37:00
% Updated: 2009.10.26 23:56:00
% ������ʾ��Ϣ

zero_cluster = find(sum(U, 2)==0);
U(zero_cluster,:) = [];

[k, data_n0] = size(U0);
[cluster_n, data_n] = size(U);
if data_n0 ~= data_n || cluster_n~=k
    fprintf('computeEntropy: �趨�����%d,����õ������%d\n', k, cluster_n);
    Entropy = 0;
    return;
end

if k <=1 || cluster_n <=1
    Entropy = 0;
    return;
end

% if (size(U0,1)~=size(U,1) || size(U0,2)~=size(U,2))
%     fprintf('computeEntropy: �������С��ͬ\n');
%     Entropy = 0;
%     return;
% end


nl = sum(U, 2);
for j = 1: k
    for l = 1: k
        njl(j, l) = U(l,:)*U0(j, :)';
    end
end

temp = zeros(k, k);
for j = 1: k
    for l = 1: k
        if (njl(j,l) ~= 0)
            temp(j,l) = njl(j,l)*log(njl(j,l)/nl(l));
        end
    end
end

Entropy = -sum(sum(temp))/(data_n*log(k));