function [DBI] = computeDBI(IDX,data,view_n)
    DBI_V = inf*ones(view_n,1);
    M.label=IDX;
    M.clusterSize=length(unique(IDX));
    if max(IDX)~=length(unique(IDX)) 
        xuzhao=unique(IDX);
        for zzxb=1:size(IDX,1)
            for xiugai=1:length(xuzhao)
                if IDX(zzxb)==xuzhao(xiugai)
                    IDX(zzxb)=xiugai;
                end
            end
        end
    end
    M.label=IDX;
    for k = 1:view_n
        M.Samples = data{k};
        DBI_V(k)=DaviesBouldin(M); 
    end
    
    tmp = 1;
    for k = 1:view_n
        tmp = tmp*DBI_V(k);
    end
    DBI = tmp^(1/view_n);
end