function [ obj ] = value( U,k,alfa,sigmal,zamda,cluster_n,DU,WEI,D)
%VALUE 此处显示有关此函数的摘要
%   此处显示详细说明
fcn=0;
for i=1:k
    for j=1:cluster_n

                 fcn = fcn+sum(sum(D(:,j),3) .*U(:,j).^2 ) +sum(sum(zamda*WEI.*abs(DU{i}-DU{i}'),2));
                  obj_w(i)=sigmal*alfa(i)*log(alfa(i));
    end 
end

obj=  fcn+sum(obj_w);
end

