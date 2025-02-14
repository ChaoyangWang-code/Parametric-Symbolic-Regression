function [theta,y_pre] = get_f_lr(p,k)
%GET_F_LR  get f1,f2,...fn by linear regression

    numGenes = numel(p.evalstr2);
    Const_pair_now=p.Const_pair_now;
    [numData,~] =size(p.ytrain);
    delta_K=p.data_i(:,1);
    Kmax=p.data_i(:,4);
    xtrain=[delta_K./k(1),Kmax./k(2)];                 
    geneOutputs = ones(numData,numGenes+2);
    check_index=1;
    for i = 1:numGenes
        ind = i + 1;
        %gene_temp=eval_exe(evalstr{i},xtrain,Const_pair_now);
        try
            % 
            gene_temp=eval([p.evalstr2{i} ';']);
        catch
            % 
            disp('error.');
        end                
        geneOutputs(:,ind)=lg(gene_temp);
        if  any(~isfinite(geneOutputs(:,ind))) || any(~isreal(geneOutputs(:,ind)))
            %disp('Error: The expression resulted in a complex number or infinity after substituting the data.');
            check_index=0;
            break
        end
    end
    
    geneOutputs(:,numGenes+2)=lg(delta_K);          
    geneOutputs_temp=geneOutputs;
    if check_index==1
        if p.ridge_on
            %% 
            geneOutputs(:,1)=[];                                                 
            theta = ridge(p.ytrain,geneOutputs,p.ridge_k,0);     
        else
            goptrans = geneOutputs';
            prj = goptrans * geneOutputs;
            theta = pinv(prj) * goptrans * p.ytrain;                     
        end        
    else
        theta=zeros(numGenes+2,1);
    end
    y_pre=geneOutputs_temp*theta;

end

