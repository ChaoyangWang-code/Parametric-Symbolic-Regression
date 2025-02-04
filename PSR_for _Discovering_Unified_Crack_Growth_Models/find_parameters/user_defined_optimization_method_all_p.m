function [loss_out,delta_kth_out,kc_out,theta_out] = user_defined_optimization_method_all_p(p)

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
    loss_MSE(1)=loss_MSE_0;


    iteration_extent=p.iteration_extent;
    iteration_limit_up=iteration_extent(1);
    iteration_limit_low=iteration_extent(2);           
    e_limit=1e-10;
    loss_limit=1e-15;
    iteration_num=0;
    e=100;
    z=ones(1,num_p);
    while ((e>e_limit)||(iteration_num<iteration_limit_low))&&(iteration_num<iteration_limit_up )&&(loss_MSE(iteration_num+1)>loss_limit)
         iteration_num=iteration_num+1;
         %% AdaDelta%%%%%%%%%%%%%%%%%%%%%%%%
        epsilon=1e-4;
        beta=0.5;
        s_g=zeros(num_p,1);
        s_t=zeros(num_p,1);
        delta_z=zeros(num_p,1);
        G=zeros(num_p,1);
        for i=1:num_p
            G(i)=1/numData.*sum( (eval(p.eq)-p.ytrain).*eval(p.diff_omega{i}).*omega0(i) );                
            s_g(i)=beta*s_g(i)+(1-beta)*G(i)^2;
            s_t(i)=beta*s_t(i)+(1-beta)*delta_z(i)^2;
            delta_z(i)=G(i)*((s_t(i)+epsilon)/(s_g(i)+epsilon))^0.5;
            %s_t(i)=beta*s_t(i)+(1-beta)*delta_z(i)^2;
        end
        z=z-delta_z';
        delta1(iteration_num)=abs(delta_z(1));
        delta2(iteration_num)=abs(delta_z(2)); 
         %% Adam%%%%%%%%%%%%%%%%%%%%%%%%
%                 downrate=0.001;
%                 epsilon=1e-7;
%                 beta=0.9;                 
%                 gamma = 0.9;
%                 s=zeros(num_p,1);
%                 m = zeros(num_p,1);
%                 G=zeros(num_p,1);
%                 for i=1:num_p
%                     G(i)=1/numData.*sum( (eval(p.eq)-p.ytrain).*eval(p.diff_omega{i}).*omega0(i) );                
%                     m(i)=gamma*m(i)+(1-gamma)*G(i);
%                     s(i)=beta*s(i)+(1-beta)*G(i)^2;
%                 end
%                 m_bar=m./(1-gamma^iteration_num);
%                 s_bar=s./(1-beta^iteration_num);
%                 delta_z=downrate./(s_bar+epsilon).^0.5.*m_bar;       
%                 z=z-delta_z';
%                 delta1(iteration_num)=abs(delta_z(1));
%                 delta2(iteration_num)=abs(delta_z(2)); 

        omega=omega0'.*z';
        f=omega(1:num_f,1);
        k=omega(num_f+1:num_p,1);        
        loss_MSE(iteration_num+1)=1/2/numData*sum((eval(p.eq)-p.ytrain).^2);
        e=abs(loss_MSE(iteration_num)-loss_MSE(iteration_num+1))/loss_MSE(iteration_num);
    end
    loss_out=loss_MSE(iteration_num+1);
    delta_kth_out=k(1);
    kc_out=k(2);
    theta_out=f';
end

