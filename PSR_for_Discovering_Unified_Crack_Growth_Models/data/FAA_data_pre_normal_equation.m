clc
clear all;
load('FAA_NASGRO_experiment_data.mat')
load('NASGRO_FIT.mat')
load('WALKER_FIT.mat')
load('FAA_data_SPLINE_FIT_240418.mat')


num_compare=size(FAA_NASGRO_experiment_data,2);
num_all_R=size(FAA_NASGRO_experiment_data_all_R,2);


for i=1:num_compare
    now_index=1;
    clear material_i_all_R
    for j=1:num_all_R
        if FAA_NASGRO_experiment_data_all_R{2,j}==i
            material_i_all_R{1,now_index}=FAA_NASGRO_experiment_data_all_R{1,j};
            now_index=now_index+1;
        end
    end
    FAA_NASGRO_experiment_data_material_R{1,i}=material_i_all_R;
end


%% 
dadN_NASGRO=cell(1,num_all_R);
for i=1:num_all_R
    data=FAA_NASGRO_experiment_data_all_R{1,i};
    data_sort=sortrows(data,1);
    delta_K=data_sort(:,1);
    dadN_experiment=data_sort(:,2);    
    R=data_sort(1,3);
    Kmax=delta_K/(1-R);
    
    SmaxS0=NASGRO_FIT{1,i}(1,1);
    alpha=NASGRO_FIT{1,i}(1,2);  
    KIc=NASGRO_FIT{1,i}(1,3);
    YS=NASGRO_FIT{1,i}(1,4);
    Ak=NASGRO_FIT{1,i}(1,5);
    Bk=NASGRO_FIT{1,i}(1,6);
    Kc=NASGRO_FIT{1,i}(1,7);
    t=NASGRO_FIT{1,i}(1,8);
    DK1=NASGRO_FIT{1,i}(1,9);
    CthP=NASGRO_FIT{1,i}(1,10);
    CthN=NASGRO_FIT{1,i}(1,11);
    p=NASGRO_FIT{1,i}(1,12);
    q=NASGRO_FIT{1,i}(1,13);
    C=NASGRO_FIT{1,i}(1,14);
    n=NASGRO_FIT{1,i}(1,15);
    %
    A0 = (0.825-0.34*alpha+0.05*alpha*alpha)*(  cos(0.5*pi*SmaxS0)  )^(1.0/alpha);
    A1 = (0.415-0.07*alpha)*SmaxS0;
    A3 = 2.0*A0 + A1 - 1.0;
    A2 = 1.0 - A0 - A1 - A3;    
    if (R>=0.0)   
        f = A0+A1*R+A2*R*R+A3*R*R*R;
	    if (f<R) 
           f=R;
        end
    elseif(R>=-2.0 && R <0.0) 
        f = A0 + A1*R;
    else
        f = A0 - 2*A1;
    end


    %delta_K1_star=DK1*(a/(a+a0))^0.5;
    delta_K1_star=DK1;

    if R>=0
        delta_Kth=delta_K1_star*((1-R)/(1-f))^(1+R*CthP)/(1-A0)^((1-R)*CthP);
    else
        delta_Kth=delta_K1_star*((1-R)/(1-f))^(1+R*CthN)/(1-A0)^(CthP-R*CthN);
    end

    %t0=2.5*(KIc/YS)^2;
    dadN_temp=C.*((1-f)/(1-R).*delta_K).^n.*(1-delta_Kth./delta_K).^p./(1-Kmax./Kc).^q;
    N_dadN=size(dadN_temp,1);
    dadN_r = [];                       
    delta_K_r=[];
    ori_position_r=[];
    dadN_i = [];
    delta_K_i=[];
    ori_position_i=[];    
    for j=1:N_dadN
        dadN_check=dadN_temp(j,1);
        if (isreal(dadN_check))&&(dadN_check>0)
            dadN_r=[dadN_r;dadN_check];
            delta_K_r=[delta_K_r;delta_K(j,1)];
            ori_position_r=[ori_position_r;j];
        else
            dadN_i=[dadN_i;dadN_check];
            delta_K_i=[delta_K_i;delta_K(j,1)];
            ori_position_i=[ori_position_i;j];
        end
    end


    
    figure(2)
    loglog(delta_K_r,dadN_r,'b-')
    hold on;grid on;
    loglog(delta_K,dadN_experiment,'r*')


    dadN_NASGRO{1,i}=[delta_K_r,dadN_r,ori_position_r];
    dadN_NASGRO{2,i}=[delta_K_i,dadN_i,ori_position_i];


end

clearvars -except dadN_NASGRO


















