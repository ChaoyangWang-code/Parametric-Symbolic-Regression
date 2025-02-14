function gp=gpdemo_crack_growth_config(gp)


rng(0, 'twister');
%run control parameters
gp.runcontrol.pop_size = 6000;                     
gp.runcontrol.num_gen = 35;		
gp.runcontrol.stage1 = 25;
gp.runcontrol.stage2 = 25;
%selection
gp.selection.tournament.size = 100;


gp.runcontrol.verbose = 1;    			%the generation frequency with which results are printed to CLI
gp.runcontrol.usecache = false;         %fitness caching: used when copying individuals in a gen



%termination
gp.fitness.plot_every_gen=0;
gp.fitness.savepath='workdir';     %example:    Windows:'C:\Users\Desktop\test250106\'   linux:'/gs/home/gp_test/test1/';
gp.fitness.terminate = true;
gp.fitness.terminate_value = 1e-20;
%
gp.fitness.loss_mode=2.5;    

gp.fitness.diff_model=2.1;
%% 
gp.fitness.loss_MSE_limit=0.03;
gp.fitness.loss_VAR_limit=1.0;
%% 
gp.fitness.auto_renew_limit=false;
gp.fitness.loss_MSE_limit_index=1.2;
gp.fitness.loss_VAR_limit_index=1.2;
gp.fitness.const_choose=[1];

gp.fitness.multidata_loss_mode=6;     
gp.fitness.cal_k_var=0;

gp.fitness.lambda=[4.0,1.5,0.1,0.5];       %
gp.fitness.entendactor=[0.9,1.1];       
gp.fitness.iteration=[5000, 20];       %
%options_qn=optimoptions("fminunc","Algorithm","quasi-newton","MaxFunctionEvaluations",2000,"OptimalityTolerance",1e-10,'Display','off');
options_qn = optimoptions('fminunc',...          
    'Algorithm', 'quasi-newton',...            
    'OptimalityTolerance', 1e-6, ...           
    'StepTolerance', 1e-10,  ...               
    'FunctionTolerance', 1e-6, ...             
    'MaxIterations', 400, ...                  
    'MaxFunctionEvaluations', 1000, ...        
    'Display', 'off');            
options_tr = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'MaxIterations',2000,'FunctionTolerance',1e-10,'MaxFunctionEvaluations',2000,'OptimalityTolerance',1e-10,'StepTolerance',1e-10,'Display','off');
options_3 = optimoptions("fmincon", "Algorithm", "interior-point", "MaxFunctionEvaluations", 2000, "OptimalityTolerance", 1e-4,'FunctionTolerance', 1e-4, 'Display', 'off');
%
gp.fitness.option1=options_qn;
gp.fitness.option2=options_tr;
gp.fitness.option3=options_3;
gp.fitness.noise_on=0;
gp.fitness.noise_level=0.8e-1; 
gp.fitness.ridge_on=1;
gp.fitness.ridge_k=0.4;

gp.fitness.theta_end_limit=[1,10];
gp.fitness.k1_k2_ini=[1,200];
gp.fitness.k1_k2_num=10;

%input configuration 
gp.nodes.inputs.num_inp = 2; 		         
gp.nodes.inputs.max_number_of_ERC = 2; 		
%quartic example doesn't need constants
gp.nodes.const.p_ERC = 0.5;		              
gp.selection.tournament.p_pareto = 1;
%gp.selection.tournament.size = gp.runcontrol.pop_size;         

%define function nodes
gp.nodes.functions.name = {'times','minus','plus','f_r'};
%genes
gp.genes.max_genes = 6;                   



%maximum depth of trees 
gp.treedef.max_depth = 6; 
 	              
%maximum depth of sub-trees created by mutation operator
gp.treedef.max_mutate_depth = 6;

%fitness function
gp.fitness.fitfun = @regressmulti_fitfun_multidata;     

%parallel
gp.runcontrol.parallel.auto = false;
gp.runcontrol.parallel.enable=false;

gp.operators.mutation.p_mutate = 0.24;    
gp.operators.crossover.p_cross = 0.64;    
gp.operators.directrepro.p_direct = 0.12; 

%load in the raw x and y data
%load demo0330data 
gp.userdata.xtrain = zeros(100,gp.nodes.inputs.num_inp); %training set (inputs)
gp.userdata.ytrain = zeros(100,1); %training set (output)
gp.userdata.xtest = zeros(100,gp.nodes.inputs.num_inp); %testing set (inputs)
gp.userdata.ytest = zeros(100,1); %testing set (output)
gp.userdata.name = 'Cherkassky function';
% load M300_preprocessing.mat
% gp.multidata=M300_preprocessing;
%load('TC4_EBM_R01_preprocessing.mat')
%gp.multidata=TC4_EBM_R01_preprocessing;
%% 
%test_data_make_NASGRO_HS;
%test_data_make_Walker;
%test_data_make_Forman;
%test_data_make_NASA;
% load('FAA_data_240417.mat')
% data_operation;
%load('FAA_data_240605.mat')
load('FAA_data4_13_combine_240628.mat')
%load('NASGRO_HS_test_no_noise_240623.mat')
data_select=false;
selection_num=1;
if data_select
    multidata{1,1}=testdata{:,selection_num};
    multidata{2,1}=1;
    R_NUM=size(testdata_temp,2);
    count=1;
    for j=1:R_NUM
        if testdata_temp{2,j}==selection_num
            all_data_without_collection{1,count}=testdata_temp{1,j};
            all_data_without_collection{2,count}=1;
            count=count+1;
        end
    end
else
    multidata=testdata;
    all_data_without_collection=testdata_temp;
end


gp.multidata=multidata;
gp.all_data_without_collection=all_data_without_collection;

