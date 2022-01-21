function [ U,alfa ] = wangdy(data,options)
tolerance=1e-3; 
MaxIter1=100;

sigmal=options.sigmal;
zamda=options.zamda;
cluster_n=options.cluster_n;
Iter=1;
k =size(data,2);
for i=1:k
[dim{i},N{i}]=size(data{i});
end
for i=1:k
[~,P0{i}] = kmeans(data{i},cluster_n);
end
for j=1:k
    for i=1:cluster_n
        dist{j}(i,:)=sqrt(sum(((data{j}-ones(dim{j}, 1)*P0{j}(i,:)).^2)'));
    end
    U=1./(dist{j}.^2.*(ones(cluster_n,1)*sum(dist{j}.^(-2))));
    U(dist{j}==0)=1;
    U=U';
end
for i=1:k
Dist{i}=dist{i}';
end
rou=2;
for i=1:k
    W_quanzhong(i)=rand(1,1);
end
W_new=W_quanzhong/sum(W_quanzhong);
alfa=W_new;
ki=6; 
for i=1:k
wei{i}=  make_affinity(data{i}, ki);
end

WEI=zeros(dim{1},dim{1});

for i =1:k
   WEI= alfa(i)*wei{i}+WEI;
end
    
% make the affinity matrix based on KNN
Nn=sum(WEI~=0,2);
CH=(WEI~=0);
CH=CH.*1;
epst=0.0001*CH ;
[dj,dk,w]=find(WEI);
au=zeros(dim{1},cluster_n);
bu=zeros(dim{1},cluster_n);

DK=cell(cluster_n,1);
for j=1:cluster_n
    DK{j}=WEI;
end
VK= DK;
pp=1:dim{1};
 for p=1:cluster_n
     pu=sparse(pp,pp,U(:,p));
     DU{p}=pu*CH;
     sst=rou* abs(DU{p} + 1/rou * DK{p}-DU{p}' - 1/rou * DK{p}' +epst);
     [~,~,vst]=find(sst);
       vst=max( 1- zamda*w./vst,0.5); 
       sst=sparse(dj,dk,vst);  
       sst1=sparse(dj,dk,1-vst);
       VK{p}= sst.*(DU{p} + 1/rou*DK{p}) +sst1.*(DU{p}' + 1/rou*DK{p}');
       DK{p}= DK{p} + rou*(DU{p}-VK{p}); 
 end 
while(Iter<MaxIter1)
    
%     for i=1:k
%         V{i}=(U.^2)'*data{i}(:,:)./((ones(size(data{i}(:,:),2),1)*sum(U.^2))'); % new center P update
%     end
%     for i=1:k
%         dist{i}=(distfcm(V{i},data{i}).^(2))';
%     end
for i=1:k
 P{i}=(U.^2)'*data{i}./(ones(N{i},1)*sum(U.^2))';
end
for j=1:k
 if size(P0{j}, 2) > 1
        for i=1:cluster_n
            Dist{j}(:,i)=sqrt(sum(((data{j}-ones(dim{1}, 1)*P0{j}(i,:)).^2)'))';         
        end
  else
        for i=1:C
            Dist{j}(:,i) = abs(P0{j}(i)-data{i});
        end
 end    
end




    for i=1:k
        D(:,:,i)=alfa(i)*Dist{i}.^2;
    end
     for i=1:cluster_n
       au(:,i)=sum(D(:,i),3) + rou/2 * Nn;  
       bu(:,i)= sum( DK{i} -rou*VK{i},2);
     end

  pu=(1+sum(bu./(2*au),2))./sum(0.5 ./ au , 2); 
  for i=1:cluster_n
      U(:,i)=(pu-bu(:,i))./(2* au(:,i)) ; 
  end
  
for p=1:cluster_n
    pu=sparse(pp,pp,U(:,i));
    DU{i}=pu*CH;
    sst=rou* abs(DU{p} + 1/rou * DK{p}-DU{p}' - 1/rou * DK{p}' +epst);
    [~,~,vst]=find(sst);
    vst=max( 1- zamda*w./vst,0.5); 
    sst=sparse(dj,dk,vst);  
    sst1=sparse(dj,dk,1-vst);
    VK{p}= sst.*(DU{p} + 1/rou*DK{p}) +sst1.*(DU{p}' + 1/rou*DK{p}');
    DK{p}= DK{p} + rou*(DU{p}-VK{p}); 
end 
     
     pp=1:dim{1};  
     for j=1:k
         for i=1:cluster_n
                 delta (j) =(sum(U(:,i).^2 .*Dist{j}(:,i).^2) +sum(sum(zamda*WEI.*abs(VK{i}-VK{i}'),2)))/sigmal;
        temp_w(j) = exp(-delta(j));

         end
     end
       
    if(sum(temp_w)==0)
        W_new = (1/k)*ones(1,k);
    else
        for i=1:k
            W_new(i)=temp_w(i)./(sum(temp_w)); 
        end
    end
alfa=W_new;
WEI=zeros(dim{1},dim{1});
for i =1:k
   WEI= alfa(i)*wei{i}+WEI;
end
    
Nn=sum(WEI~=0,2);
CH=(WEI~=0);
CH=CH.*1;
epst=0.0001*CH ;
[dj,dk,w]=find(WEI);

      cost(Iter)=value( U,k,alfa,sigmal,zamda,cluster_n,DU,WEI,D);
    if Iter>1
        if abs(cost(Iter)-cost(Iter-1))<tolerance
          fprintf('µü´ú´ÎÊý:  %d\n',Iter);
            break;
        end
    end 
    P0=P;
       Iter=Iter+1;  
end
end


    
  