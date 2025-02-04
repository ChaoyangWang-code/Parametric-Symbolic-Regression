function [loss_multidata] = cal_loss_multidata(gp,fitness_best_constC_all_data,numGenes,num_const,lambda)
    theta=fitness_best_constC_all_data(:,4+num_const:4+num_const+numGenes+1);           %theta: numGenes+2
    lambda1=lambda(1);
    lambda2=lambda(2);
    lambda3=lambda(3);
    lambda4=lambda(4);

    multidata_loss=gp.fitness.multidata_loss_mode;

    if multidata_loss==1
%% 
        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=inf;                                   
        else
            loss_multidata=max(fitness_best_constC_all_data(:,1));                 
        end


    elseif multidata_loss==2
%%         
        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=inf;
        else
            Ndata=size(theta,1);
            m=mean(theta);
            ss=std(theta);
            max_of_var=max(ss);
            if max_of_var>lambda1
                loss_multidata=inf;
            else
                loss_of_fit=lambda2*max(fitness_best_constC_all_data(:,1));
                loss_of_mean=lambda3/Ndata*(sum(m.^2));
                loss_of_var=lambda4*sum(ss);
                loss_multidata=loss_of_fit+loss_of_mean+loss_of_var;
            end        
    %         loss_of_mean=lambda1/Ndata*(sum(m.^2));
    %         loss_of_var=lambda2*(max(v));
    %         loss=loss_of_fit+loss_of_mean+loss_of_var;
           %fitness_at_constC_all_data(:,1)=loss;
        end
    
    elseif multidata_loss==3
        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=inf;
        else
            Ndata=size(theta,1);
            m=mean(theta);                       
            ss=std(theta);                            
            max_of_var=max(ss);
            loss_of_fit=max(fitness_best_constC_all_data(:,1));
            %loss_of_var=sum(ss);
            loss_of_var=max_of_var;

            if gp.fitness.gen_count_now<gp.runcontrol.stage1                     
                loss_multidata=loss_of_fit;             
            else

                if loss_of_fit>gp.fitness.loss_MSE_limit          
                    loss_multidata=inf;                                     
                else
                    loss_multidata=loss_of_var;
                end
            end
        end

    elseif multidata_loss==4
        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=inf;
        else
            Ndata=size(theta,1);                               
            loss_of_fit=max(fitness_best_constC_all_data(:,1));
            k_of_AIC=size(theta,2);

            if gp.fitness.gen_count_now<gp.runcontrol.stage1                      
                loss_multidata=loss_of_fit;             
            else

                if loss_of_fit>gp.fitness.loss_MSE_limit          
                    loss_multidata=inf;                                     
                else
                    for i=1:Ndata
                        n_of_AIC=size(gp.multidata{1,i},1);      
                        AIC(i)=2*k_of_AIC+n_of_AIC*log(fitness_best_constC_all_data(i,1));
                        AICc(i)=AIC(i)+2*k_of_AIC*(k_of_AIC+1)/(n_of_AIC-k_of_AIC-1);
                    end                    
                    loss_multidata=max(AICc);
                    if loss_multidata<0
                        loss_multidata=-1/loss_multidata;   
                    else
                        loss_multidata=inf;
                    end

                end

            end

        end
    elseif multidata_loss==5  

        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=[inf,inf,inf];
        else
            Ndata=size(theta,1);
            m=mean(theta);                        
            ss=std(theta);                           
            max_of_var=max(ss);
            loss_of_fit=max(fitness_best_constC_all_data(:,1));
            %loss_of_var=sum(ss);
            loss_of_var=max_of_var;
            loss_of_parameter_num=numGenes+2;

            %loss_of_fit
            if gp.fitness.gen_count_now<=gp.runcontrol.stage1                       
                loss_multidata=[loss_of_fit,loss_of_var,loss_of_parameter_num];             
            end
            %loss_of_var，and loss_of_fit<gp.fitness.loss_MSE_limit
            if (gp.fitness.gen_count_now>gp.runcontrol.stage1)&&(gp.fitness.gen_count_now<=gp.runcontrol.stage2)   
                if loss_of_fit<gp.fitness.loss_MSE_limit          
                    %loss_multidata=loss_of_var;
                    loss_multidata=[loss_of_fit,loss_of_var,loss_of_parameter_num];
                else
                    loss_multidata=[inf,inf,inf];      
                    %loss_multidata=inf;                                     
                end   
            end
            %loss_of_parameter_num，and loss_of_fit<gp.fitness.loss_MSE_limit，loss_of_var<gp.fitness.loss_VAR_limit
            if gp.fitness.gen_count_now>gp.runcontrol.stage2   
                if (loss_of_fit<gp.fitness.loss_MSE_limit)&&(loss_of_var<gp.fitness.loss_VAR_limit)   
                    loss_multidata=[loss_of_fit,loss_of_var,loss_of_parameter_num];
                    %loss_multidata=loss_of_parameter_num;
                else
                    loss_multidata=[inf,inf,inf];
                    %loss_multidata=inf;                                     
                end   

            end
        end
%% 
    elseif multidata_loss==6    
        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=[inf,inf,inf];
        else
            Ndata=size(theta,1);
            m=mean(theta);                       
            %ss=std(theta);                           
            ss=var(theta);                            
            if gp.fitness.cal_k_var
                k1=fitness_best_constC_all_data(:,2);
                k2=fitness_best_constC_all_data(:,3);
                %max_of_var=max([ss,std(k1),std(k2)]);
                max_of_var=max([ss,var(k1),var(k2)]);
            else
                max_of_var=max(ss);
            end
            loss_of_fit=max(fitness_best_constC_all_data(:,1));

            if Ndata==1
                loss_of_var=0.1;
            else
                %loss_of_var=sum(ss);
                loss_of_var=max_of_var;
            end

            loss_of_parameter_num=numGenes+2;
            loss_multidata=[loss_of_fit,loss_of_var,loss_of_parameter_num];
        end
    else

        error('wrong multidata_loss')

    end



end


