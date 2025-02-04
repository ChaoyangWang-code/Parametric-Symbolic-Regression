function [fitness_out,gp,theta,ypredtrain,fitnessTest,ypredtest,pvals,r2train,r2test,r2val,geneOutputs,geneOutputsTest,geneOutputsVal]=regressmulti_fitfun_multidata(evalstr_in,gp)
%REGRESSMULTI_FITFUN Fitness function for multigene symbolic regression.
%
%   This is the default fitness function for multigene symbolic regression
%   in GPTIPS.
%
%   [FITNESS,GP] = REGRESSMULTI_FITFUN(EVALSTR,GP) returns the FITNESS of
%   the symbolic expression(s) in the cell array EVALSTR using information
%   contained in the GP data struct. Here, FITNESS is the root mean squared
%   prediction error (RMSE) on the training data set.
%
%   [FITNESS,GP,THETA,YPREDTRAIN,FITNESS_TEST,YPREDTEST,PVALS,R2TRAIN,R2TEST,R2VAL]
%   = REGRESSMULTI_FITFUN(EVALSTR,GP) may be used post-run to return the
%   gene coefficients THETA, the prediction of the model on the training
%   data YPREDTRAIN, the RMSE fitness value FITNESS_TEST on the test data
%   set, the prediction of the model on the test data YPREDTEST, the
%   statistical p-values for bias and model terms are returned as PVALS
%   (PVALS only computed if the Statistics Toolbox is present, otherwise an
%   empty variable is returned). Additionally, coefficients of
%   determination (R^2) are returned as R2TRAIN, R2TEST and R2VAL.
%
%   Remarks:
%
%   Each observation of the response variable y is assumed to be an unknown
%   non-linear function of the corresponding observations of the predictor
%   variables x1,..xn.
%
%   Training data:
%
%   The GPTIPS configuration file should populate the following required
%   fields for the training data assuming 'Ntrain' observations on the
%   input and output data. GP.USERDATA.XTRAIN should be a (Ntrain X n)
%   matrix where the ith column contains the Ntrain observations of the ith
%   input variable xi. GP.USERDATA.YTRAIN should be a (Ntrain x 1) vector
%   containing the corresponding observations of the response variable y.
%
%   Testing data:
%
%   The following fields are optional and may be used, post-run, to see how
%   well evolved models generalise to an unseen test data set with Ntest
%   observations. They do not affect the model building process.
%   GP.USERDATA.XTEST should be a (Ntest X n) matrix where the ith column
%   contains the Ntest observations of the ith input variable xi.
%   GP.USERDATA.YTEST should be a (Ntest x 1) vector containing the
%   corresponding observations of the response variable y.
%
%   How multigene symbolic regression works:
%
%   In multigene symbolic regression, each prediction of y is formed by the
%   weighted output of each of the trees/genes in the multigene individual
%   plus a bias term. The number (M) and structure of the trees is evolved
%   automatically during a GPTIPS run (subject to user defined
%   constraints).
%
%   i.e. ypredtrain = c0 + c1*tree1 + ... + cM*treeM
%
%   where c0 = bias term
%         c1,..,cM are the weights
%         M is the number of genes/trees comprising the current individual
%
%   The weights (i.e. regression coefficients) are automatically determined
%   by a least squares procedure for each multigene individual and are
%   stored in GP.FITNESS.RETURNVALUES for future use.
%
%   Remarks:
%
%   Because the GP structure is modified within this function (i.e. the
%   field GP.FITNESS.RETURNVALUES is used to store the computed weighting
%   coefficients for each gene) the GP structure must be returned as an
%   output argument.
%
%   This fitness function is used for multigene symbolic regression for
%   GPDEMO2, GPDEMO3 and GPDEMO4 (the configuration files for these are
%   GPDEMO2_CONFIG.M and GPDEMO3_CONFIG.M respectively) but it can and
%   should be used for the user's own non-linear regression problems.
%
%   Copyright (c) 2009-2015 Dominic Searson
%   Copyright (c) 2023-2025 Chaoyang Wang
%   GPTIPS 2
%
%   See also REGRESSMULTI_FITFUN_VALIDATE, GPDEMO2_CONFIG, GPDEMO3_CONFIG,
%   GPDEMO4_CONFIG, GPDEMO2, GPDEMO3

%defaults in case of early exit
theta=[];ypredtrain=[];fitnessTest=[];ypredtest=[];
r2train=[];r2test=[];r2val=[];geneOutputs=[];geneOutputsTest=[];
geneOutputsVal=[];
%% for parallel
multidata=gp.multidata;
bootSample=gp.userdata.bootSample;
bootSampleSize=gp.userdata.bootSampleSize;
run_completed=gp.state.run_completed;
force_compute_theta=gp.state.force_compute_theta;
iteration_extent=gp.fitness.iteration;
%% 
numGenes = numel(evalstr_in); 
num_const=0;
for i=1:numGenes
    open_sq_br = strfind(evalstr_in{i},'c');
    num_const = num_const+numel(open_sq_br);
end
const_choose=gp.fitness.const_choose;
size_const_choose=size(const_choose,2);
const_all_situations=size_const_choose^(num_const);
const_all=zeros(const_all_situations,num_const);
for i=1:const_all_situations
    for j=1:num_const
        if j==1
            index=mod(i-1,size_const_choose)+1;      
            const_all(i,num_const+1-j)=const_choose(index);
        else
            index=floor(mod(i-1,size_const_choose^j)/(size_const_choose^(j-1)))+1;
            const_all(i,num_const+1-j)=const_choose(index);
        end
    end
end

Ndata=size(multidata,2);
% gp.fitness.num_const=num_const;
% gp.fitness.numGenes=numGenes;
% gp.fitness.Ndata=Ndata;
%% 
loss_mode=gp.fitness.loss_mode;
diff_model=gp.fitness.diff_model;
if loss_mode>=2
    [eq,diff_omega,contains_k1,contains_k2]=diff_F(evalstr_in,num_const,diff_model);

    pata = 'k(\d+)';
    patb = 'c(\d+)';
    patc = 'f(\d+)';
    diff_omega = regexprep(diff_omega,pata,'k($1)');
    diff_omega = regexprep(diff_omega,patb,'Const_pair_now($1)');
    diff_omega = regexprep(diff_omega,patc,'f($1)');
    eq = regexprep(eq,pata,'k($1)');
    eq = regexprep(eq,patb,'Const_pair_now($1)');
    eq = regexprep(eq,patc,'f($1)');

    parameter_gp.diff_omega=diff_omega;
    parameter_gp.eq=eq;
    parameter_gp.contains_k1=contains_k1;
    parameter_gp.contains_k2=contains_k2;
%     parameter_gp.dF_dk1=dF_dk1;
%     parameter_gp.dF_dk2=dF_dk2;
end

if contains_k1&&contains_k2
   k1k2_Correlation=k1k2_Correlation_test(evalstr_in);           
   parameter_gp.k1k2_Correlation=k1k2_Correlation;
end

%% 
if loss_mode>=10  
    if diff_model==2.1
        DELTA_y_DELTA_k1=cell(1,Ndata);
        DELTA_y_DELTA_k2=cell(1,Ndata);
        R=cell(1,Ndata);
        k1=cell(1,Ndata);
        k2=cell(1,Ndata);
        for data_j=1:Ndata
            D_i=multidata{1,data_j};
            [num_of_data_i,~]=size(D_i(:,1));
            R{data_j}=D_i(1,3);
            DELTA_y_DELTA_k1{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   D_i(2:num_of_data_i,1)-D_i(1:num_of_data_i-1,1)   );
            DELTA_y_DELTA_k2{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   D_i(2:num_of_data_i,4)-D_i(1:num_of_data_i-1,4)   );
            k1{data_j}=0.5*(  D_i(1:num_of_data_i-1,1) +D_i(2:num_of_data_i,1)  );
            k2{data_j}=0.5*(  D_i(1:num_of_data_i-1,4) +D_i(2:num_of_data_i,4)  );
        end

    elseif diff_model==2.2
        DELTA_y_DELTA_z1=cell(1,Ndata);
        DELTA_y_DELTA_z2=cell(1,Ndata);
        R=cell(1,Ndata);
        z1=cell(1,Ndata);
        z2=cell(1,Ndata);
        for data_j=1:Ndata
            D_i=multidata{1,data_j};
            [num_of_data_i,~]=size(D_i(:,1));
            R{data_j}=D_i(1,3);
            DELTA_y_DELTA_z1{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   lg(D_i(2:num_of_data_i,1))-lg(D_i(1:num_of_data_i-1,1))   );
            DELTA_y_DELTA_z2{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   lg(D_i(2:num_of_data_i,4))-lg(D_i(1:num_of_data_i-1,4))   );
            z1{data_j}=0.5*(  lg(D_i(1:num_of_data_i-1,1)) +lg(D_i(2:num_of_data_i,1))  )  ;
            z2{data_j}=0.5*(  lg(D_i(1:num_of_data_i-1,4)) +lg(D_i(2:num_of_data_i,4))  )  ;
        end

    end
end

%% 
if contains_k2
    KC=linspace(log10(gp.fitness.k1_k2_ini(1)),log10(gp.fitness.k1_k2_ini(2)),gp.fitness.k1_k2_num);  
else
    KC=log10(1);    
end
if contains_k1
    DELTA_KTH=linspace(log10(gp.fitness.k1_k2_ini(1)),log10(gp.fitness.k1_k2_ini(2)),gp.fitness.k1_k2_num); 
else
    DELTA_KTH=log10(1);
end
if contains_k1&&contains_k2
    if k1k2_Correlation
        DELTA_KTH=log10(1);
    end    
end
                       
DELTA_KTH=10.^DELTA_KTH;
KC=10.^KC;


%% 
for const_parameter_i=1:const_all_situations       
    check_at_const_parameter_i=1;                
    Const_pair_now=const_all(const_parameter_i,:);
    parameter_gp.Const_pair_now=Const_pair_now;
    fitness_best_constC_all_data=zeros(Ndata,numGenes+6+num_const)*inf;
%%     
    for data_i=1:Ndata                                             
        if  check_at_const_parameter_i==0
            break;
        end
        if loss_mode>=10         
            if diff_model==2.1
                parameter_gp.k1=k1{data_i};
                parameter_gp.k2=k2{data_i};
                parameter_gp.DELTA_y_DELTA_k1=DELTA_y_DELTA_k1{data_i};
                parameter_gp.DELTA_y_DELTA_k2=DELTA_y_DELTA_k2{data_i};   
            elseif diff_model==2.2
                parameter_gp.z1=z1{data_i};
                parameter_gp.z2=z2{data_i};
                parameter_gp.DELTA_y_DELTA_z1=DELTA_y_DELTA_z1{data_i};
                parameter_gp.DELTA_y_DELTA_z2=DELTA_y_DELTA_z2{data_i};   
            end  
        end
        %parameter_gp.R=R{data_i};
        parameter_gp.data_i=multidata{1,data_i};
        [num_of_data_i,~]=size(multidata{1,data_i}(:,1));
        parameter_gp.num_of_data_i=num_of_data_i;        
        %%
        %log_delta_K_i_min=min(lg(multidata{1,data_i}(:,1)));                     
        %log_delta_K_i_max=max(lg(multidata{1,data_i}(:,4)));
%         DELTA_KTH=(log_delta_K_i_min-0.4):0.05:(log_delta_K_i_min+0.1);
%         KC=(log_delta_K_i_max-0.2):0.05:(log_delta_K_i_max+0.3);
        %dis_1=(log_delta_K_i_min-log10(1))/5;
        %dis_2=(log10(200)-log_delta_K_i_max)/5;
        % DELTA_KTH=(log10(1)):dis_1:(log_delta_K_i_min-0.001);                 
        % KC=(log_delta_K_i_max+0.001):dis_2:(log10(200));                         
        % DELTA_KTH=linspace(log10(gp.fitness.k1_k2_ini(1)),log10(gp.fitness.k1_k2_ini(2)),gp.fitness.k1_k2_num);               
        % KC=linspace(log10(gp.fitness.k1_k2_ini(1)),log10(gp.fitness.k1_k2_ini(2)),gp.fitness.k1_k2_num);          
        % DELTA_KTH=10.^DELTA_KTH;
        % KC=10.^KC;
        %DELTA_KTH=flip(DELTA_KTH);
        %KC=flip(KC);                
        %DELTA_KTH=[1.82,DELTA_KTH];
        %KC=[KC,142.5];    
        [~,n1]=size(DELTA_KTH);
        [~,n2]=size(KC);
        N_material=n1*n2;
        Parameters_meterial=zeros(N_material,2);
        for i=1:N_material
            index2=mod(i-1,n2)+1; 
            index1=floor(mod(i-1,n1*n2)/(n2))+1;
            Parameters_meterial(i,1)= DELTA_KTH(index1);   
            Parameters_meterial(i,2)= KC(index2); 
        end
        loss_at_Kth_KC_all=zeros(N_material,1)*inf;
        fitness_at_Kth_KC_all=cell(N_material,1);
        geneOutputs=cell(N_material,1);       
        %% 
        parameter_gp.evalstr1=evalstr_in;
        parameter_gp.fit_option1=gp.fitness.option1;
        parameter_gp.fit_option2=gp.fitness.option2;
        parameter_gp.fit_option3=gp.fitness.option3;
        parameter_gp.noise_on=gp.fitness.noise_on;
        parameter_gp.noise_level=gp.fitness.noise_level; 
        parameter_gp.ridge_on=gp.fitness.ridge_on;        
        parameter_gp.ridge_k=gp.fitness.ridge_k;
        parameter_gp.theta_end_limit=gp.fitness.theta_end_limit;
        p_temp=parameter_gp;
        evalstr=evalstr_in;
        evalstr_test=evalstr_in;     
        %%      
        for Kth_KC_i=1:N_material                  
        %parfor Kth_KC_i=1:N_material
            theta=zeros(numGenes+2,1);
            fitness_temp_at_Kth_KC_i=zeros(1,5+num_const+numGenes);
            check_fitness=1;
            parameter_K=Parameters_meterial(Kth_KC_i,:);                        
            delta_kth=parameter_K(1);
            kc=parameter_K(2);
            xtrain=[];
            ytrain=[];
%             %%%delta_Kth/delta_K, delta_K/delta_Kth, Kmax/Kc, Kc/Kmax, R, log10(delta_K)
%             xtrain=[delta_kth./multidata{1,data_i}(:,1),multidata{1,data_i}(:,1)./delta_kth,multidata{1,data_i}(:,4)./kc,kc./multidata{1,data_i}(:,4),multidata{1,data_i}(:,3)];
            %delta_K/delta_Kth, Kmax/Kc, log10(delta_K)
            xtrain=[multidata{1,data_i}(:,1)./delta_kth,multidata{1,data_i}(:,4)./kc];
            ytrain=log10(multidata{1,data_i}(:,2));
            % process evalstr with regex to allow direct access to data matrices
            pat1 = 'x(\d+)';
            pat2 = 'c(\d+)';
            evalstr = regexprep(evalstr,pat2,'Const_pair_now($1)');
            evalstr = regexprep(evalstr,pat1,'xtrain(:,$1)');
            p_temp.evalstr2=evalstr;
            y = ytrain;
            [numData,~] =size(ytrain);
            %set up a matrix to store the tree outputs plus a bias column of ones
            geneOutputs{Kth_KC_i,1} = ones(numData,numGenes+2);
            %% 
            for i = 1:numGenes
                ind = i + 1;
                %gene_temp=eval_exe(evalstr{i},xtrain,Const_pair_now);
                try
                    % 
                    gene_temp=eval([evalstr{i} ';']);
                catch
                    % 
                    disp('An error occurred.');
                end

                geneOutputs{Kth_KC_i,1}(:,ind)=lg(gene_temp);
                %eval(['geneOutputs{Kth_KC_i,1}(:,ind)=' evalstr{i} ';']);
                %check for nonsensical answers and break out early with an 'inf' if so
                if  any(~isfinite(geneOutputs{Kth_KC_i,1}(:,ind))) || any(~isreal(geneOutputs{Kth_KC_i,1}(:,ind)))
                    %fitness = Inf;
                    %gp.fitness.returnvalues{gp.state.current_individual} = [];
                    %return
                    fitness_temp_at_Kth_KC_i=[Inf,delta_kth,kc,Const_pair_now,zeros(1,numGenes+2),0];  
                    check_fitness=0;
                    break
                end
            end
            if check_fitness==1
                loss_cal_optimize_in=1;
                geneOutputs{Kth_KC_i,1}(:,numGenes+2)=lg(multidata{1,data_i}(:,1));                
                    %only calc. weighting coeffs during an actual run or if forced
                if ~run_completed || force_compute_theta
            
                    %set gp.userdata.bootSample to true to resample data for weights computation
                    %prepare LS matrix
                    if bootSample
                        sampleInds = bootsample(geneOutputs{Kth_KC_i,1},bootSampleSize);
                        goptrans = geneOutputs{Kth_KC_i,1}(sampleInds,:)';
                        prj = goptrans * geneOutputs{Kth_KC_i,1}(sampleInds,:);
                        ysample = y(sampleInds);
                    else
                        goptrans = geneOutputs{Kth_KC_i,1}';
                        prj = goptrans * geneOutputs{Kth_KC_i,1};
                    end
            
                    %calculate tree weight coeffs using SVD based least squares
                    %normal equation
                    try
                        if bootSample
                            theta = pinv(prj) * goptrans * ysample;
                            if p_temp.ridge_on
                                %% 
                                geneOutputs{Kth_KC_i,1}(:,1)=[];
                                theta = ridge(ysample,geneOutputs{Kth_KC_i,1},p_temp.ridge_k,0);             
                            end                                 
                        else
                            theta = pinv(prj) * goptrans * y;
                            if p_temp.ridge_on
                                %% 
                                geneOutputs{Kth_KC_i,1}(:,1)=[];
                                theta = ridge(y,geneOutputs{Kth_KC_i,1},p_temp.ridge_k,0);             
                            end                                  
                        end                        
                    catch
                        %theta = [];
                        %fitness = Inf;
                        %gp.fitness.returnvalues{gp.state.current_individual} = [];
                        %return;
                         loss_cal_optimize_in=0;
                         fitness_temp_at_Kth_KC_i=[Inf,delta_kth,kc,Const_pair_now,zeros(1,numGenes+2),0];  
                    end
            
                    %assign bad fitness if any coeffs NaN or Inf
                    if any(isinf(theta)) || any(isnan(theta))
                        %theta = [];
                        %fitness = Inf;
                        %gp.fitness.returnvalues{gp.state.current_individual} = [];
                        %return;
                         loss_cal_optimize_in=0;
                        fitness_temp_at_Kth_KC_i=[Inf,delta_kth,kc,Const_pair_now,zeros(1,numGenes+2),0];                  
                    end
            
                    %write coeffs to returnvalues field for storage
                    %gp.fitness.returnvalues{gp.state.current_individual} = theta;
                else %if post-run, get stored coeffs from return value field
                    %theta = gp.fitness.returnvalues{gp.state.current_individual};
                    error('no theta!');
                end    
                if loss_cal_optimize_in==1  
       
                    p_temp.theta=theta;
                    p_temp.parameter_K=parameter_K;
                    p_temp.iteration_extent=iteration_extent;
                    p_temp.geneOutputs=geneOutputs{Kth_KC_i,1};
                    p_temp.ytrain=ytrain;
                    %loss=loss_cal(loss_mode,p_temp);
                    [loss,delta_kth,kc,theta,opti_mode_used]=loss_cal_optimize(loss_mode,p_temp,DELTA_KTH,KC);
                    %% 
                    materials_N=size(gp.multidata,2);                  
                    data_N=size(gp.all_data_without_collection,2);    
                    data_all=gp.all_data_without_collection;
                    data_all_collection=gp.multidata;                
                    material_i_Kmax_min=min(data_all_collection{1,data_i}(:,4));
                    material_i_Kmax_max=max(data_all_collection{1,data_i}(:,4));
                    material_i_R_min=min(data_all_collection{1,data_i}(:,3));
                    %material_i_R_max=max(data_all_collection{1,data_i}(:,3));    
                    material_i_R_max=0.99;
                    if material_i_R_min~=material_i_R_max
                        material_i_Kmax_pre = linspace(material_i_Kmax_min, material_i_Kmax_max,100);
                        %material_i_delta_K_pre = linspace(material_i_delta_K_min, material_i_delta_K_max,100);
                        material_i_R_pre = linspace(material_i_R_min, material_i_R_max,100);  
                        [material_i_Kmax_PREDICTION, material_i_R_PREDICTION] = meshgrid(material_i_Kmax_pre, material_i_R_pre);
                        % xtest{1}=material_i_delta_K_PREDICTION./delta_kth;
                        % xtest{2}= material_i_delta_K_PREDICTION./(1-material_i_R_PREDICTION)./kc;
                        material_i_delta_K_PREDICTION=material_i_Kmax_PREDICTION.*(1-material_i_R_PREDICTION);
                        xtest{1}=material_i_delta_K_PREDICTION./delta_kth;
                        xtest{2}= material_i_Kmax_PREDICTION./kc;                        
                        num_of_predata_i=size(material_i_Kmax_pre,2);
                        pat11 = 'x(\d+)';
                        pat22 = 'c(\d+)';
                        evalstr_test = regexprep(evalstr_test,pat22,'Const_pair_now($1)');
                        evalstr_test = regexprep(evalstr_test,pat11,'xtest{$1}'); 
    
                        geneOutputs_test = theta(1).*ones(num_of_predata_i,num_of_predata_i);
                        for j_test = 1:numGenes
                            gene_temp=eval([evalstr_test{j_test} ';']);
                            geneOutputs_test=geneOutputs_test+theta(j_test+1).*lg(gene_temp);
                        end 
                        geneOutputs_test=geneOutputs_test+theta(2+numGenes).*lg(material_i_delta_K_PREDICTION);       
                        if isreal(geneOutputs_test) && ~any(isinf(geneOutputs_test),'all')
                            test_index = 1;     
                            fitness_temp_at_Kth_KC_i=[loss,delta_kth,kc,Const_pair_now,theta',opti_mode_used];     %3+num_const+numGenes+2+1
                        else
                            test_index = 0; 
                            fitness_temp_at_Kth_KC_i=[inf,delta_kth,kc,Const_pair_now,theta',opti_mode_used];     %3+num_const+numGenes+2+1
                        end
    
                    else 
                        fitness_temp_at_Kth_KC_i=[loss,delta_kth,kc,Const_pair_now,theta',opti_mode_used];     %3+num_const+numGenes+2+1
                    end

                end
 
                %fitness_temp_at_Kth_KC_i=[loss,delta_kth,kc,Const_pair_now,theta',opti_mode_used];     %3+num_const+numGenes+2+1
            end
            loss_at_Kth_KC_all(Kth_KC_i,1)=fitness_temp_at_Kth_KC_i(1);
            fitness_at_Kth_KC_all{Kth_KC_i,1}=fitness_temp_at_Kth_KC_i;    
            if fitness_temp_at_Kth_KC_i(1)<gp.fitness.terminate_value
                break
            end

        end  
        [loss_min_at_Kth_KC_all,pvals]=min(loss_at_Kth_KC_all);      
        if  isinf(loss_min_at_Kth_KC_all)
            check_at_const_parameter_i=0;          
        end
        if check_at_const_parameter_i==1
            fitness_best_constC_all_data(data_i,1:numGenes+6+num_const)=fitness_at_Kth_KC_all{pvals,1};     
        end
    end 
    loss_multidata=cal_loss_multidata(gp,fitness_best_constC_all_data,numGenes,num_const,gp.fitness.lambda);
    fitness_all_const{1}{const_parameter_i,1}=fitness_best_constC_all_data;
    fitness_all_const{2}(const_parameter_i,1)=loss_multidata(1);   
    fitness_all_const{3}(const_parameter_i,1)=loss_multidata(2);   
    fitness_all_const{4}(const_parameter_i,1)=loss_multidata(3);   
end
if gp.fitness.gen_count_now<=gp.runcontrol.stage1        
    [~,pvals]=min(fitness_all_const{1,2});       
end
if (gp.fitness.gen_count_now>gp.runcontrol.stage1)&&(gp.fitness.gen_count_now<=gp.runcontrol.stage2)   
    [~,pvals]=min(fitness_all_const{1,3}); 
end
if gp.fitness.gen_count_now>gp.runcontrol.stage2 
    [~,pvals]=min(fitness_all_const{1,4});       
end

fitness_out=[fitness_all_const{2}(pvals,1),fitness_all_const{3}(pvals,1),fitness_all_const{4}(pvals,1)];

fitness_temp=fitness_all_const{1,1}{pvals,1};
if loss_mode>=2
   fitness_return={fitness_temp,numGenes,num_const,eq,evalstr_in};
elseif loss_mode==1
   fitness_return={fitness_temp,numGenes,num_const,evalstr_in};
end
gp.fitness.returnvalues = fitness_return;
%foldername=['W:\test231012' sprintf('%d',inpChoice)];
