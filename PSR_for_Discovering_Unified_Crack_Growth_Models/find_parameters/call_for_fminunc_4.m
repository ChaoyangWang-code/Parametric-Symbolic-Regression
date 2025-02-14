function [fit,G] = call_for_fminunc_4(z)

global global_p_temp

numGenes = numel(global_p_temp.evalstr2);
Const_pair_now=global_p_temp.Const_pair_now;
f0=global_p_temp.theta';
k0=global_p_temp.parameter_K;
omega0=[f0,k0];                       
%omega=omega0';
[numData,~] =size(global_p_temp.ytrain);
delta_K=global_p_temp.data_i(:,1);
Kmax=global_p_temp.data_i(:,4);
eq=global_p_temp.eq;
ytrain=global_p_temp.ytrain;
num_p=size(global_p_temp.diff_omega,2);
num_f=num_p-2;

%% update k,f
omega=omega0'.*z';
f=omega(1:num_f,1);
k=omega(num_f+1:num_p,1);        


%% Calculate fit, G
fit=1/2/numData*sum((eval(eq)-ytrain).^2);

if nargout > 1 % gradient required
    for i=1:num_p
        G(i,1)=1/numData.*sum( (eval(eq)-ytrain).*eval(global_p_temp.diff_omega{i}).*omega0(i) );                
    end
end


end

