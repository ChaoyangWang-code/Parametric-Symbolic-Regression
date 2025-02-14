function [fit,G] = call_for_fminunc_5(z)
    
    global global_p_temp
    global search_index
    global loss_MSE  
    global loss_min            
    global z_now
    global theta_now
    global z_min
    global theta_min

    numGenes = numel(global_p_temp.evalstr2);
    Const_pair_now=global_p_temp.Const_pair_now;
    f=global_p_temp.theta';
    k0=global_p_temp.parameter_K;
    %omega0=[f,k];                       
    %omega=omega0';
    [numData,~] =size(global_p_temp.ytrain);
    delta_K=global_p_temp.data_i(:,1);
    Kmax=global_p_temp.data_i(:,4);
    eq=global_p_temp.eq;
    ytrain=global_p_temp.ytrain;
    num_p=numGenes+4;
    num_f=num_p-2;
    %% update k
    k=z.*k0;
    if global_p_temp.contains_k2==0
        k(2)=k0(2);          
    end
    if global_p_temp.contains_k1==0
        k(1)=k0(1);          
    end
    if global_p_temp.contains_k1&&global_p_temp.contains_k2
        if global_p_temp.k1k2_Correlation
            k(1)=k0(1);          
        end    
    end

    
    [theta,y_pre]=get_f_lr(global_p_temp,k);
    %% update f
    f=theta';
    fit=1/2/numData*sum((y_pre-ytrain).^2);
    
    if nargout > 1 % gradient required
        for i=1:2
            G(i,1)=1/numData.*sum( (y_pre-ytrain).*eval(global_p_temp.diff_omega{num_f+i}).*k0(i) );               
        end
    end
    search_index=search_index+1;
    loss_MSE(search_index)=fit;
    z_now=z;
    theta_now=theta;
    if fit<loss_min
        loss_min=fit;
        z_min=z;
        theta_min=theta;
    end


end

