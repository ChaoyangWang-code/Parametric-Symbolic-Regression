function [loss_out,delta_kth_out,kc_out,theta_out] = user_defined_optimization_method_k1k2(p)

    numGenes = numel(p.evalstr2);
    Const_pair_now=p.Const_pair_now;
    f=p.theta';
    k=p.parameter_K;
    k0=k;
    [numData,~] =size(p.ytrain);
    delta_K=p.data_i(:,1);
    Kmax=p.data_i(:,4);
    num_p=size(p.diff_omega,2);
    num_f=num_p-2;
    loss_MSE_0=sum((eval(p.eq)-p.ytrain).^2)/2/numData;    
    

    f_now=f;
    k_now=k;

    loss_MSE(1)=loss_MSE_0;
    loss_MSE_min=loss_MSE_0;
    iteration_extent=p.iteration_extent; 
    e_limit=1e-6;
    iteration_limit_up=iteration_extent(1);
    iteration_limit_low=iteration_extent(2);
    loss_limit=1e-15;
    iteration_num=0;
    e=100;
    z=[1,1];

    while ((e>e_limit)||(iteration_num<iteration_limit_low))&&(iteration_num<iteration_limit_up )&&(loss_MSE(iteration_num+1)>loss_limit)

        %% 
%                 downrate=1.0;
%                 m=0;
%                 gama=0.9;
%                 G=zeros(2,1);
%                 for i=1:2
%                     G(i)=1/numData.*sum( (eval(p.eq)-p.ytrain).*eval(p.diff_omega{num_f+i}).*k0(i) );               
%                 end                
%                 m=gama*m+downrate.*G;
%                 z=z-m';
%                 delta1(iteration_num+1)=abs(m(1));
%                 delta2(iteration_num+1)=abs(m(2));

        %% Nesterov%%%%%%%%%%%%%%%%%%%%%%%%
%                 downrate=1;
%                 m=0;
%                 gama=0.9;
%                 G=zeros(2,1);
%                 for i=1:2
%                     G(i)=1/numData.*sum( (eval(p.eq)-p.ytrain).*eval(p.diff_omega{num_f+i}).*k0(i) );                
%                 end      
%                 m=gama*m+downrate.*G.*(z'-gama*m);  
%                 z=z-m';
%                 delta1(iteration_num+1)=abs(m(1));
%                 delta2(iteration_num+1)=abs(m(2));

        %% RMSProp%%%%%%%%%%%%%%%%%%%%%%%%        
%                 epsilon=1e-6;
%                 downrate=0.001;  
%                 beta=0.99;
%                 m=0;
%                 gama=0.9;                        
%                 s=zeros(2,1);
%                 G=zeros(2,1);       
%                 for i=1:2
%                     G(i)=1/numData.*sum( (eval(p.eq)-p.ytrain).*eval(p.diff_omega{num_f+i}).*k0(i) );                
%                     s(i)=beta*s(i)+(1-beta)*G(i)^2;
%                 end                
%                 delta_z=(downrate./(s+epsilon).^0.5.*G)';
%                 z=z-delta_z;
%                 delta1(iteration_num+1)=abs(delta_z(1));
%                 delta2(iteration_num+1)=abs(delta_z(2));

        %% AdaDelta%%%%%%%%%%%%%%%%%%%%%%%%
        epsilon=1e-4;
        beta=0.9;
        s_g=zeros(2,1);
        s_t=zeros(2,1);
        delta_z=zeros(2,1);
        G=zeros(2,1);
        for i=1:2
            G(i)=1/numData.*sum( (eval(p.eq)-p.ytrain).*eval(p.diff_omega{num_f+i}).*k0(i) );                
            s_g(i)=beta*s_g(i)+(1-beta)*G(i)^2;
            s_t(i)=beta*s_t(i)+(1-beta)*delta_z(i)^2;
            delta_z(i)=G(i)*((s_t(i)+epsilon)/(s_g(i)+epsilon))^0.5;
            %s_t(i)=beta*s_t(i)+(1-beta)*delta_z(i)^2;
        end
        z=z-delta_z';
        delta1(iteration_num+1)=abs(delta_z(1));
        delta2(iteration_num+1)=abs(delta_z(2)); 
        %% Adam%%%%%%%%%%%%%%%%%%%%%%%%
%                 downrate=0.1;
%                 epsilon=1e-8;
%                 beta=0.8;         
%                 gamma = 0.8;
%                 s=zeros(2,1);
%                 m = zeros(2,1);
%                 G=zeros(2,1);
%                 for i=1:2
%                     G(i)=1/numData.*sum( (eval(p.eq)-p.ytrain).*eval(p.diff_omega{num_f+i}).*k0(i) );                
%                     m(i)=gamma*m(i)+(1-gamma)*G(i);
%                     s(i)=beta*s(i)+(1-beta)*G(i)^2;
%                 end
%                 m_bar=m./(1-gamma^(iteration_num+1));
%                 s_bar=s./(1-beta^(iteration_num+1));
%                 delta_z=downrate./(s_bar+epsilon).^0.5.*m_bar;       
%                 z=z-delta_z';
%                 delta1(iteration_num+1)=abs(delta_z(1));
%                 delta2(iteration_num+1)=abs(delta_z(2)); 
        %% Nadam%%%%%%%%%%%%%%%%%%%%%%%%
%                 downrate=0.2;
%                 epsilon=1e-6;
%                 beta=0.9;         
%                 gamma = 0.9;
%                 s=zeros(2,1);
%                 m = zeros(2,1);
%                 G=zeros(2,1);
%                 for i=1:2
%                     G(i)=1/numData.*sum( (eval(p.eq)-p.ytrain).*eval(p.diff_omega{num_f+i}).*k0(i) );               
%                     m(i)=gamma*m(i)+(1-gamma)*G(i);
%                     s(i)=beta*s(i)+(1-beta)*G(i)^2;
%                 end
%                 G_bar=G./(1-gamma^(iteration_num+1));
%                 m_bar=m./(1-gamma^(iteration_num+1));
%                 s_bar=s./(1-beta^(iteration_num+1));
%                 m_mean=(1-gamma)*G_bar+gamma*m_bar;
%                 delta_z=downrate./(s_bar+epsilon).^0.5.*m_mean;       
%                 z=z-delta_z';
%                 delta1(iteration_num+1)=abs(delta_z(1));
%                 delta2(iteration_num+1)=abs(delta_z(2)); 

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         k=k0.*z;

        %% update loss
        xtrain=[delta_K./k(1),Kmax./k(2)];
        geneOutputs = ones(numData,numGenes+2);
        for i = 1:numGenes
            ind = i + 1;
            %gene_temp=eval_exe(evalstr{i},xtrain,Const_pair_now);
            try
                % 
                gene_temp=eval([p.evalstr2{i} ';']);
            catch
                % 
                disp('An error occurred.');
            end

            geneOutputs(:,ind)=lg(gene_temp);
            if  any(~isfinite(geneOutputs(:,ind))) || any(~isreal(geneOutputs(:,ind)))
                %fitness = Inf;
                %gp.fitness.returnvalues{gp.state.current_individual} = [];
                %return
                loss=inf;
                break
            end
        end
        geneOutputs(:,numGenes+2)=lg(delta_K); 

%                 goptrans = geneOutputs';
%                 prj = goptrans * geneOutputs;
%                 theta = pinv(prj) * goptrans * p.ytrain;

        geneOutputs(:,1)=[];
        kk=0;
        theta = ridge(p.ytrain,geneOutputs,kk,0);
        f=theta';
        iteration_num=iteration_num+1;
        loss_MSE(iteration_num+1)=1/2/numData*sum((eval(p.eq)-p.ytrain).^2);
        e=abs(loss_MSE(iteration_num)-loss_MSE(iteration_num+1))/loss_MSE(iteration_num);
        if loss_MSE(iteration_num+1)<loss_MSE_min
            loss_MSE_min=loss_MSE(iteration_num+1);
            f_now=f;
            k_now=k;                    
        end

    end

    loss_out=loss_MSE_min;
    delta_kth_out=k_now(1);
    kc_out=k_now(2);
    theta_out=f_now';


end

