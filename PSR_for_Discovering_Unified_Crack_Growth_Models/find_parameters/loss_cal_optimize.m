function [loss_out,delta_kth_out,kc_out,theta_out,opti_mode_used] = loss_cal_optimize(loss_mode,p,DELTA_KTH,KC)

    if loss_mode==1
        %% 
        [numData,~] =size(p.ytrain);
        ypredtrain = p.geneOutputs * p.theta;                          
        err = p.ytrain - ypredtrain;
        loss_out=sqrt(((err'*err)/numData));
        theta_out=p.theta';
        delta_kth_out=p.parameter_K(1);
        kc_out=p.parameter_K(2);

    elseif loss_mode>=2   
        %% 

        if loss_mode==2.1

            [loss_out,delta_kth_out,kc_out,theta_out] = user_defined_optimization_method_all_p(p);
       
        elseif loss_mode==2.2

            [loss_out,delta_kth_out,kc_out,theta_out] = user_defined_optimization_method_k1k2(p);

        elseif loss_mode==2.4
            %% 
            numGenes = numel(p.evalstr2);
            Const_pair_now=p.Const_pair_now;
            f=p.theta';
            k=p.parameter_K;
            k0=k;
            f0=f;
            omega0=[f0,k0];                     
            [numData,~] =size(p.ytrain);
            delta_K=p.data_i(:,1);
            Kmax=p.data_i(:,4);
            num_p=size(p.diff_omega,2);
            num_f=num_p-2;
            loss_MSE_0=sum((eval(p.eq)-p.ytrain).^2)/2/numData;   
            global loss_MSE
            loss_MSE(1)=loss_MSE_0;

            global global_p_temp
            global_p_temp=p;
            z0=ones(1,num_p);
            fun = @rosenbrockwithgrad_4;
            try
                % 
                [z,loss_now,exitflag,output] = fminunc(fun,z0,p.fit_option1);
            catch
                % 
                try
                    % 
                    [z,loss_now,exitflag,output] = fminunc(fun,z0,p.fit_option2);
                catch
                    % 
                    loss_now=loss_MSE(1);
                    z=z0;
                end
            end
            clear global_p_temp  
            clear loss_MSE

            if loss_now<loss_MSE(1)
                loss_out_temp=loss_now;
                omega=omega0'.*z';
                f=omega(1:num_f,1);
                k=omega(num_f+1:num_p,1);      

                delta_kth_out=k(1);
                kc_out=k(2);
                theta_out=f;        
            else
                loss_out_temp=loss_MSE(1);
                delta_kth_out=k0(1);
                kc_out=k0(2);
                theta_out=f0';               
            end
            %% theta_end check
            if (theta_out(end)>theta_end_limit(1)) &&(theta_out(end)<theta_end_limit(2))
                loss_out=loss_out_temp;
            else
                loss_out=inf;                     
            end

        elseif loss_mode==2.5
            numGenes = numel(p.evalstr2);
            Const_pair_now=p.Const_pair_now;
            f=p.theta';
            k=p.parameter_K;
            k0=k;
            f0=f;
            omega0=[f0,k0];                 
            [numData,~] =size(p.ytrain);
            delta_K=p.data_i(:,1);
            Kmax=p.data_i(:,4);
            num_p=size(p.diff_omega,2);
            num_f=num_p-2;
            loss_MSE_0=sum((eval(p.eq)-p.ytrain).^2)/2/numData;   
            z0 = [1 1];

            global loss_MSE
            loss_MSE=[];
            loss_MSE(1)=loss_MSE_0;      
            global loss_min
            loss_min=[];
            loss_min=loss_MSE_0; 

            global z_now         
            z_now=[];
            z_now=z0;         
            global theta_now
            theta_now=[];
            theta_now=p.theta;  
            global z_min           
            z_min=[];
            z_min=z0;         
            global theta_min
            theta_min=[];
            theta_min=p.theta;  

            global search_index
            search_index=1;
            global global_p_temp
            global_p_temp=p;



            loss_now=loss_MSE_0;
            fun = @rosenbrockwithgrad_5;
            try
                % 
                [~,~,exitflag,output] = fminunc(fun,z0,p.fit_option1);
                opti_mode_used=1;
            catch
                % 
                try
                    % 
                    loss_MSE=[];
                    loss_MSE(1)=loss_MSE_0;     
                    z_now=[];
                    z_now=z0;  
                    theta_now=[];
                    theta_now=p.theta;     
                    loss_min=[];
                    loss_min=loss_MSE_0;         
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    theta_min=p.theta;  
                    [~,~,exitflag,output] = fminunc(fun,z0,p.fit_option2);       
                    opti_mode_used=2;
                catch
                    % 
                    %loss_min=loss_MSE(end);
                    %z=z0;
                    loss_MSE=[];
                    loss_MSE(1)=loss_MSE_0;     
                    z_now=[];
                    z_now=z0;  
                    theta_now=[];
                    theta_now=p.theta;     
                    loss_min=[];
                    loss_min=loss_MSE_0;         
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    theta_min=p.theta;                      
                    opti_mode_used=0;
                end
            end

            loss_out_temp=loss_min;
            k=z_min.*k0;
            f=theta_min';
            delta_kth_out=k(1);
            kc_out=k(2);            
            theta_out=f';

            clear global_p_temp
            clear search_index
            clear loss_MSE  
            clear loss_min            
            clear z_now
            clear theta_now
            clear z_min
            clear theta_min
            %% 
            %
            pass_index=1;
            if theta_out(end)<p.theta_end_limit(1) 
                pass_index=0;
            elseif theta_out(end)>p.theta_end_limit(2)
                pass_index=0;
            elseif (delta_kth_out<0)||(kc_out<0)
                pass_index=0;
            % elseif delta_kth_out>kc_out
            %     pass_index=0;                        
            end

            if pass_index==1
                loss_out=loss_out_temp;
            else
                loss_out=inf;
            end

        elseif loss_mode==2.6
            %% 
            numGenes = numel(p.evalstr2);
            Const_pair_now=p.Const_pair_now;
            f=p.theta';
            k=p.parameter_K;
            k0=k;
            f0=f;
            omega0=[f0,k0];                       
            [numData,~] =size(p.ytrain);
            delta_K=p.data_i(:,1);
            Kmax=p.data_i(:,4);
            num_p=size(p.diff_omega,2);
            num_f=num_p-2;
            loss_MSE_0=sum((eval(p.eq)-p.ytrain).^2)/2/numData;   
            z0 = [1 1];
            lb = [min(DELTA_KTH), min(KC)];  
            ub = [max(DELTA_KTH), max(KC)];  


            global loss_MSE
            loss_MSE=[];
            loss_MSE(1)=loss_MSE_0;      
            global loss_min
            loss_min=[];
            loss_min=loss_MSE_0; 

            global z_now           
            z_now=[];
            z_now=z0;         
            global theta_now
            theta_now=[];
            theta_now=p.theta;  
            global z_min           
            z_min=[];
            z_min=z0;         
            global theta_min
            theta_min=[];
            theta_min=p.theta;  

            global search_index
            search_index=1;
            global global_p_temp
            global_p_temp=p;

            loss_now=loss_MSE_0;
            fun = @rosenbrockwithgrad_5;
            try
                % 
                [~,~,exitflag,output]= fmincon(fun, z0, [], [], [], [], lb, ub, [], p.fit_option3);       
                opti_mode_used=3;
            catch
                % 
                try
                    % 
                    loss_MSE=[];
                    loss_MSE(1)=loss_MSE_0;     
                    z_now=[];
                    z_now=z0;  
                    theta_now=[];
                    theta_now=p.theta;     
                    loss_min=[];
                    loss_min=loss_MSE_0;         
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    theta_min=p.theta;  
                    [~,~,exitflag,output] = fminunc(fun,z0,p.fit_option1);       
                    opti_mode_used=1;
                catch
                    % 
                    %loss_min=loss_MSE(end);
                    %z=z0;
                    loss_MSE=[];
                    loss_MSE(1)=loss_MSE_0;     
                    z_now=[];
                    z_now=z0;  
                    theta_now=[];
                    theta_now=p.theta;     
                    loss_min=[];
                    loss_min=loss_MSE_0;         
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    theta_min=p.theta;                      
                    opti_mode_used=0;
                end
            end

            loss_out_temp=loss_min;
            k=z_min.*k0;
            f=theta_min';
            delta_kth_out=k(1);
            kc_out=k(2);            
            theta_out=f';

            clear global_p_temp
            clear search_index
            clear loss_MSE  
            clear loss_min            
            clear z_now
            clear theta_now
            clear z_min
            clear theta_min
            %% 
            %
            pass_index=1;
            % if theta_out(end)<p.theta_end_limit(1) 
            %     pass_index=0;
            % elseif theta_out(end)>p.theta_end_limit(2)
            %     pass_index=0;
            if (delta_kth_out<0)||(kc_out<0)
                pass_index=0;
            % elseif delta_kth_out>kc_out
            %     pass_index=0;                        
            end

            if pass_index==1
                loss_out=loss_out_temp;
            else
                loss_out=inf;
            end

        end
   

    else
        %%
        error('wrong loss_mode!')

    end
end

%load('B:\Desktop\matlab.mat')
%loss_cal_optimize(loss_mode,p)