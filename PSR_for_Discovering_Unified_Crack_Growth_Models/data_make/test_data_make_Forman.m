%% make NASGRO(H-S)data

%% 
data_all=[
-0.2	2.79	2.12	79
0.3	    2.79	2.12	79
0.5	    2.79	2.12	79
0.2	    3.12	2.54	62
0.4	    3.12	2.54	62
0.6	    3.12	2.54	62
0.8	    3.12	2.54	62
];


noise_level =gp.fitness.noise_level;
%figure
for i=1:7
    R=data_all(i,1);
    D=data_all(i,2)*1e-10;
    p=data_all(i,3);
    Kc=data_all(i,4);


    delta_K_temp=log10(1+0.1):0.02:log10(60*(1-R));
    delta_K=10.^(delta_K_temp);
    rn=size(delta_K,2);
    
    K_max=delta_K/(1-R);
    dadN=D.*( delta_K ).^p./(Kc./K_max-1);

    if gp.fitness.noise_on==1
        noise = noise_level * randn(1, rn);       
        dadN_temp=log10(dadN)+noise;
    else
        dadN_temp=log10(dadN);
    end

    dadN=10.^(dadN_temp);
    %dadN=D.*( (delta_K-k1) ./(1-K_max/k2).^0.5  ).^p+noise;
    plot(log10(delta_K),log10(dadN),'r.')
    hold on; grid on;
    n_delta_K=size(delta_K,2);
    Rn=R*ones(n_delta_K,1);
    testdata_temp{i}=[delta_K',dadN',Rn,K_max'];
end
for i = 1:3
    testdata_temp{2,i} = 1;
end
for i = 4:7
    testdata_temp{2,i} = 2;
end



disp(['number of test data:    ' num2str(size(testdata,2))]);







