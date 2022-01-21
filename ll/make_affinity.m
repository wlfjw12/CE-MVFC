function [WEI]=  make_affinity(Ldata, k)
cdata=Ldata;
[Num,Np]=size(cdata); 
 k3=k+1;
 w3=sparse(zeros(Num,Num));
 w=sparse(zeros(Num,Num));
  for i=1:Num
  temp=cdata(i,:);
  tem=temp(ones(Num,1),:);
  w1=sqrt(sum((cdata-tem).^2,2));
  [dist position] = sort(w1);
  w3(i,position(2:k3))=1;
  w1=exp(-w1./0.5);
  w(i,:)=w3(i,:).*(w1');
  end

 WEI=w+0.001*w';
