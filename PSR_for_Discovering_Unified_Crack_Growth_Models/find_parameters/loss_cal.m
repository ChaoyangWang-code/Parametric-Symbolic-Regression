function [loss] = loss_cal(loss_mode,p)

    if loss_mode==1
        %%
        [numData,~] =size(p.ytrain);
        ypredtrain = p.geneOutputs * p.theta;                          
        err = p.ytrain - ypredtrain;
        loss=sqrt(((err'*err)/numData));
    elseif loss_mode>=2
        if loss_mode==2.1
        %%
            %%%delta_y_delta_k1=dF_dk1+dF_dk2/(1-R)
            %%%delta_y_delta_k2=dF_dk2+dF_dk1*(1-R)
            theta=p.theta;
            parameter_K=p.parameter_K;
            Const_pair_now=p.Const_pair_now;
            num_of_data_i=p.num_of_data_i;
            k1=p.k1;
            k2=p.k2;
            R=p.R;
            delta_y_delta_k1=eval(p.dF_dk1)+eval(p.dF_dk2)./(1-R);
            delta_y_delta_k2=eval(p.dF_dk2)+eval(p.dF_dk1).*(1-R);
            loss1=1/num_of_data_i*(  sum( lg( 1+abs(p.DELTA_y_DELTA_k1-delta_y_delta_k1) ) )   );
            loss2=1/num_of_data_i*(  sum( lg( 1+abs(p.DELTA_y_DELTA_k2-delta_y_delta_k2) ) )   );
            loss=max([loss1,loss2]); 

        elseif loss_mode==2.2
        %%
            %%%delta_y_delta_z1=dF_dz1+dF_dz2
            %%%delta_y_delta_z2=dF_dz2+dF_dz1
            theta=p.theta;
            parameter_K=p.parameter_K;
            Const_pair_now=p.Const_pair_now;
            num_of_data_i=p.num_of_data_i;
            z1=p.z1;
            z2=p.z2;
            R=p.R;
            delta_y_delta_z1=eval(p.dF_dk1)+eval(p.dF_dk2);                    
            delta_y_delta_z2=eval(p.dF_dk2)+eval(p.dF_dk1);
%             loss1=1/num_of_data_i*(  sum( lg( 1+abs(p.DELTA_y_DELTA_z1-delta_y_delta_z1) ) )   );
%             loss2=1/num_of_data_i*(  sum( lg( 1+abs(p.DELTA_y_DELTA_z2-delta_y_delta_z2) ) )   );
            loss1=1/num_of_data_i*(  sum( lg( 1+(p.DELTA_y_DELTA_z1-delta_y_delta_z1).^2 ) )   );
            loss2=1/num_of_data_i*(  sum( lg( 1+(p.DELTA_y_DELTA_z2-delta_y_delta_z2).^2 ) )   );
            loss=max([loss1,loss2]); 
        end

    else
        %%
        error('wrong loss_mode!')

    end
end

%loss_cal(loss_mode,p)