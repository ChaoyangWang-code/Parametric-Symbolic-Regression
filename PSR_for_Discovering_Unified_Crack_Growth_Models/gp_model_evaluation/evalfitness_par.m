function gp = evalfitness_par(gp,gen_count)
%EVALFITNESS_PAR Calls the user specified fitness function (parallel version).
%
%   GP = EVALFITNESS_PAR(GP) evaluates the the fitnesses of individuals
%   stored in the GP structure and updates various other fields of GP
%   accordingly.
%
%   This has the same functionality as EVALFITNESS but makes use of the
%   Mathworks Parallel Computing Toolbox to distribute the fitness
%   computations across multiple cores.
%
%   Copyright (c) 2009-2015 Dominic Searson
%   Copyright (c) 2023-2025 Chaoyang Wang
%   GPTIPS 2
%
%   See also TREE2EVALSTR, EVALFITNESS
gp.fitness.gen_count_now=gen_count;
popSize = gp.runcontrol.pop_size;
popSize_new=gp.state.new_pop_num;
evalstrs = cell(popSize_new,1);
complexityMeasure = gp.fitness.complexityMeasure;
complexities = zeros(popSize_new,1);
fitfun = gp.fitness.fitfun;
fitvals = zeros(popSize_new,1);
returnvals = cell(popSize_new,1);
values_alldata=cell(popSize_new,1);

if gp.runcontrol.usecache
    usecache = true;
else
    usecache = false;
end



%% #############################################
%for i = 1:gp.state.new_pop_num
parfor i = 1:gp.state.new_pop_num
    %assign copy of gp inside parfor loop to a temp struct
    tempgp = gp;
    %update state of temp variable to index of the individual that is about to
    %be evaluated
    tempgp.state.current_individual = i;
    
    if usecache && tempgp.fitness.cache.isKey(i)
        
        cache = tempgp.fitness.cache(i);
        complexities(i) = cache.complexity;
        fitvals(i) = cache.value;
        returnvals{i,1} = cache.returnvalues;
        
    else
        
        %process coded trees into evaluable matlab expressions
        evalstrs{i} = tree2evalstr(tempgp.pop{i},gp);
        
        %store complexity of individual (either number of nodes or tree
        % "expressional complexity").
        if complexityMeasure
            complexities(i) = getcomplexity(tempgp.pop{i});
        else
            complexities(i) = getnumnodes(tempgp.pop{i});
        end
        
        %evaluate gp individual using fitness function
        [fitness,tempgp] = feval(fitfun,evalstrs{i},tempgp);       
        returnvals{i,1} = tempgp.fitness.returnvalues;
        values_alldata{i,1} = fitness;

        %% loss
        %loss_fit
        if tempgp.fitness.gen_count_now<=tempgp.runcontrol.stage1
            fitvals(i)=fitness(1);       
        end
        %loss_of_var
        if (tempgp.fitness.gen_count_now>tempgp.runcontrol.stage1)&&(tempgp.fitness.gen_count_now<=tempgp.runcontrol.stage2)   
           if fitness(1)<tempgp.fitness.loss_MSE_limit
                fitvals(i)=fitness(2);
           else
                fitvals(i)=inf;
           end
        end
        %loss_of_parameter_num
        if tempgp.fitness.gen_count_now>tempgp.runcontrol.stage2 
            if (fitness(1)<tempgp.fitness.loss_MSE_limit)&&(fitness(2)<tempgp.fitness.loss_VAR_limit)
                fitvals(i)=fitness(3);
            else
                fitvals(i)=inf;
            end
            %fitvals(i)=fitness(3);
        end
        %fitvals(i)=max(fitness);
        %fitvals(i)=sum( fitness)/size(fitness,2);
    end

end %end of parfor loop
%attach returned values to original GP structure
%gp.fitness.values = max(values_all_pop);

%% 
new_results_this_gen=cell(gp.state.new_pop_num,5);   %pop，returnvalues，values_alldata，gp.fitness.values，gp.fitness.complexity
new_results_this_gen(:,1)=gp.pop(1:gp.state.new_pop_num,:);      
new_results_this_gen(:,2)=returnvals;
new_results_this_gen(:,3)=values_alldata;    %
for i = 1:gp.state.new_pop_num
    new_results_this_gen{i,4}=fitvals(i);
    new_results_this_gen{i,5}=complexities(i);
end

if gen_count==1
    gp.result_all_history=new_results_this_gen;
    clear new_results_this_gen;
else
    result_all_history_temp=[gp.result_all_history;new_results_this_gen];
    gp.result_all_history=result_all_history_temp;
    clear new_results_this_gen;
    clear result_all_history_temp;
end

for i = 1:(popSize-gp.state.new_pop_num)
    position_of_old_pop=gp.positions(i);
    index=gp.state.new_pop_num+i;
    gp.pop{index,1}=gp.result_all_history{position_of_old_pop,1};
    returnvals{index,1}=gp.result_all_history{position_of_old_pop,2};
    values_alldata{index,1} = gp.result_all_history{position_of_old_pop,3};

    if gp.fitness.gen_count_now<=gp.runcontrol.stage1
        fitvals(index,1)=values_alldata{index,1}(1);       
    end

    if (gp.fitness.gen_count_now>gp.runcontrol.stage1)&&(gp.fitness.gen_count_now<=gp.runcontrol.stage2) 
       if values_alldata{index,1}(1)<gp.fitness.loss_MSE_limit
            fitvals(index,1)=values_alldata{index,1}(2);
       else
            fitvals(index,1)=inf;
       end
        %fitvals(index,1)=values_alldata{index,1}(2);
    end

    if gp.fitness.gen_count_now>gp.runcontrol.stage2 
       if (values_alldata{index,1}(1)<gp.fitness.loss_MSE_limit)&&(values_alldata{index,1}(2)<gp.fitness.loss_VAR_limit)
            fitvals(index,1)=values_alldata{index,1}(3);
       else
            fitvals(index,1)=inf;
       end      
        %fitvals(index,1)=values_alldata{index,1}(3);
    end

    complexities(index,1)=gp.result_all_history{position_of_old_pop,5};
end


gp.fitness.values = fitvals;
gp.fitness.values_alldata=values_alldata;
gp.fitness.returnvalues = returnvals;
gp.fitness.complexity = complexities;
gp.state.current_individual = popSize;

str=num2str(gen_count);
foldername_gen=[gp.fitness.savepath, str];

if ~exist(foldername_gen, 'dir')
    mkdir(foldername_gen);
end
savePath=[gp.fitness.savepath, str,'\gen',str, '.mat'];
save(savePath);






%show_results_gen_i(gp,foldername_gen,str)
%load('B:\Desktop\matlab240619.mat')









