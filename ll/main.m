clear;
clc;

load('mul_orl.mat')
name='mul_orl';
Y=labels;
data=data';
vn = length(data);
Z=importdata('bestH_ORL.mat');
Z=Z';
data{vn+1}=Z;
cluster = max(Y); 
sigmals =2.^(-1);          
zamdas=10.^(0);            
for i=1:vn
   
    data{i}=mapminmax(data{i},0,1);
    data_DBI{i} = data{i};
end
results = [];
for sigmal= sigmals
    for zamda = zamdas 
        options.sigmal=sigmal;
        options.zamda = zamda;
        options.cluster_n= cluster;
        try
        for ll=1:10
            tic
            [ U,alfa] = wangdy(data,options);
            [~,U1_index] = max(U');
            N = size(data{1},1);
            U1 = defuzzy( U1_index,N,cluster);
            exteral_indicators = evaluate_cluster_indicators(Y,U1_index');
            dbi = computeDBI(U1_index',data_DBI,vn);
            indi_temp(ll,:)=[exteral_indicators,dbi];
            time(ll)=toc; 
        end
        catch err
    disp(err);
   warning('something error');
    break;
        end     
        indicator_mean=mean(indi_temp);
        indicator_std=std(indi_temp);
        time_mean=mean(time);
        time_std=std(time);
        entropy_mean=indicator_mean(1);
        entropy_std=indicator_std(1);
        purity_mean=indicator_mean(2);
        purity_std=indicator_std(2);
        NMI_mean=indicator_mean(3);
        NMI_std=indicator_std(3);
        RI_mean=indicator_mean(4);
        RI_std=indicator_std(4);
        recall_mean=indicator_mean(5);
        recall_std=indicator_std(5);
        precision_mean=indicator_mean(6);
        precision_std=indicator_std(6);
        f_measure_mean=indicator_mean(7);
        f_measure_std=indicator_std(7);
        fmi_mean=indicator_mean(8);
        fmi_std=indicator_std(8);
        DBI_mean=indicator_mean(9);
        DBI_std=indicator_std(9);
        results = [results;[sigmal,zamda,NMI_mean,NMI_std,RI_mean,RI_std,DBI_mean,DBI_std,precision_mean,precision_std,purity_mean,purity_std,entropy_mean,entropy_std,recall_mean,recall_std,f_measure_mean,f_measure_std,fmi_mean,fmi_std, time_mean,time_std]];
        save(strcat('result',name),'results');

    end
end
    
