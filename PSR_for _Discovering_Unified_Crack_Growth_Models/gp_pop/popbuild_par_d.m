function gp = popbuild_par_d(gp,pop_all)
%POPBUILD Build next population of individuals.
%
%   GP = POPBUILD(GP) uses the current population (stored in GP.POP) and
%   their fitnesses (stored in GP.FITNESS.VALUES) and, optionally,
%   complexities to create the next generation of individuals.
%
%   Copyright (c) 2009-2015 Dominic Searson
%   Copyright (c) 2023-2025 Chaoyang Wang
%   GPTIPS 2
%
%   See also INITBUILD

%initialise new population


tour_ind = 1:gp.runcontrol.pop_size;
tour_ind=tour_ind';
tour_fitness = gp.fitness.values;
tour_comp  = gp.fitness.complexity;
%perform fast pareto sort on tournament members
if gp.fitness.minimisation
    mo = [tour_fitness tour_comp];
else
    mo = [ (-1 *tour_fitness) tour_comp];
end
%remove any 'infs' from consideration and recall if necessary
infs= isinf(mo(:,1));
mo(infs,:)=[];
if isempty(mo)
    error('All fitness are inf');
end
tour_ind(infs,:)=[];
rank = ndfsort_rank1(mo);        
rank1_tour_ind_temp = tour_ind(rank == 1);     
%num_rank1 = numel(rank1_tour_ind_temp);
pareto_fit=tour_fitness(rank1_tour_ind_temp,1);
pareto_comp=tour_comp(rank1_tour_ind_temp,1);
%%
p_all=[pareto_fit,pareto_comp];
pu=unique(p_all,'rows','stable');              
pu = sortrows(pu, 1);          
pu_num=size(pu,1);
rank1_tour_ind=[];
for ui=1:pu_num
    row_to_find = pu(ui,:);
    matches = ismember(p_all, row_to_find, 'rows');
    [row, ~] = find(matches);
    rank1_tour_ind(ui,1)=rank1_tour_ind_temp(row(1),1);
end
num_rank1 = pu_num;
newPop_pareto=cell(num_rank1,1);
for pi=1:num_rank1
    pareto_index=rank1_tour_ind(pi);
    newPop_pareto{pi,1}=gp.pop{pareto_index,1};
end

%% elite pop
%the number of new members to be constructed after elitism is accounted for
num2build = floor( (1 - gp.selection.elite_fraction) * gp.runcontrol.pop_size);
num2skim = gp.runcontrol.pop_size - num2build;

if (num2skim+num_rank1)>gp.runcontrol.pop_size
    num2skim=gp.runcontrol.pop_size-num_rank1;
end

num2skim1=round(num2skim/2);
num2skim2=num2skim-num2skim1;

newPop_skim = cell(num2skim,1);

%skim off the existing elite individuals and stick them on the end of the
%new population. When sorting there may be many individuals with the same
%fitness. However, some of these have fewer nodes/lower complexity than
%others, so skim off the first one as the existing member with the fewest
%nodes. This should exert a degree of parsimony pressure.

%get indices of best num2skim individuals
[~,sortIndex] = sort(gp.fitness.values);

%if maximising need to flip vectors
if ~gp.fitness.minimisation
    sortIndex = flipud(sortIndex);
end


for g=1:num2skim1
    elite_of_all_pop_id=sortIndex(g);
    newPop_skim{g,1} = gp.pop{elite_of_all_pop_id,1};
    loss_elite_gen(g,1)=gp.fitness.values(elite_of_all_pop_id);
end
gp.fitness.loss_elite_gen=loss_elite_gen;
gp.fitness.pop_elite_gen=newPop_skim;


for g=1:num2skim2
    
    oldIndex = sortIndex(g);
    
    %for the first individual to skim off, pick the 'best' member with the
    %lowest complexity
    if g==1
        bestInds = find( gp.fitness.values <= gp.fitness.entendactor*(gp.fitness.values(oldIndex)) );
        [~,oldIndex] = sort(gp.fitness.complexity(bestInds));
        oldIndex = bestInds(oldIndex(1)); %if more than one with same complexity just pick the first
    end
    
    newIndex = g+num2skim1;
    copiedIndividual = gp.pop{oldIndex,1};
    
    %cache fitness, complexity etc.
    if gp.runcontrol.usecache
        cachedData.complexity = gp.fitness.complexity(oldIndex,1);
        cachedData.returnvalues = gp.fitness.returnvalues{oldIndex,1};
        cachedData.value = gp.fitness.values(oldIndex,1);
        gp.fitness.cache(newIndex) = cachedData;
    end
    
    newPop_skim{newIndex,1} = copiedIndividual;
    
end

%% 
%the number of new members to be constructed after elitism is accounted for
num2build = gp.runcontrol.pop_size-num_rank1-num2skim;
newPop_temp = cell(num2build,1);
%parameter shortcuts
p_mutate = gp.operators.mutation.p_mutate;
p_direct = gp.operators.directrepro.p_direct;
pmd = p_mutate + p_direct;

%reset cache
if gp.runcontrol.usecache
    remove(gp.fitness.cache,gp.fitness.cache.keys);
end

% num2build=100000
% newPop_temp = cell(gp.runcontrol.pop_size,1);
%update gen counter
gp.state.count = gp.state.count + 1;
directPop=zeros(num2build,1);
%% #############################################
%for buildCount=1:num2build
parfor buildCount=1:num2build
    tempgp = gp;
    check_index=0;
    while check_index==0
        %probabilistically select a genetic operator
        p_gen = rand;
        
        if p_gen < p_mutate  %select mutation

            eventType = 1;
            parentIndex = selection(tempgp);  %pick the population index of a parent individual using selection operator
            parent = tempgp.pop{parentIndex};
            [check_index,newPop_out]=popbuild_eT1(parent,tempgp);
            if check_index==1
                newPop_temp{buildCount,1}={newPop_out,[]};
                break;
            end

        elseif p_gen < pmd   %direct reproduction
            eventType = 2;
            parentIndex = selection(tempgp);  %pick a parent
            parent = tempgp.pop{parentIndex};
            %copy to new population
            newPop_temp{buildCount,1}={parent,[]};
            directPop(buildCount,1)=parentIndex;
            %store fitness etc of copied individual if cache enabled
%             if tempgp.runcontrol.usecache
%                 cachedData.complexity = tempgp.fitness.complexity(parentIndex,1);
%                 cachedData.returnvalues = tempgp.fitness.returnvalues{parentIndex,1};
%                 cachedData.value = tempgp.fitness.values(parentIndex,1);
%                 tempgp.fitness.cache(buildCount) = cachedData;
%             end
            break;
        else                 %crossover
            eventType = 3;
%             try
%                 % 
%                 [check_index,newPop_out]=popbuild_eT3(tempgp);
%             catch
%                 % 
%                 disp('An error occurred.');
%             end

            [check_index,newPop_out]=popbuild_eT3(tempgp);
            if check_index==1
                newPop_temp{buildCount,1}=newPop_out;
                break;
            end
        end
    end

end

%newPop_temp
newPop_build = cell(num2build,1);
newPop_num=0;
for i=1:num2build
    newPop_pair=newPop_temp{i,1};
    if ~isempty(newPop_temp{i,1}{1,2})
        newPop_num=newPop_num+1;
        if newPop_num>num2build
            break;
        end
        newPop_build{newPop_num,1}=newPop_temp{i,1}{1,1};
        newPop_num=newPop_num+1;
        if newPop_num>num2build
            break;
        end
        newPop_build{newPop_num,1}=newPop_temp{i,1}{1,2};
    else 

        if directPop(i)>0          
            pindex=directPop(i);
            newPop_num=newPop_num+1;
            if newPop_num>num2build
                break;
            end            
            newPop_build{newPop_num,1}=newPop_temp{i,1}{1,1};
            %store fitness etc of copied individual if cache enabled
            if gp.runcontrol.usecache
                cachedData.complexity = gp.fitness.complexity(pindex,1);
                cachedData.returnvalues = gp.fitness.returnvalues{pindex,1};
                cachedData.value = gp.fitness.values(pindex,1);
                gp.fitness.cache(newPop_num) = cachedData;
            end
        else                     
            newPop_num=newPop_num+1;
            if newPop_num>num2build
                break;
            end                   
            newPop_build{newPop_num,1}=newPop_temp{i,1}{1,1};
        end

    end
end


newPop = vertcat(newPop_pareto, newPop_skim, newPop_build);


%%%%
gp.state.pareto_in_thisgen=num_rank1;
gp.state.elite_in_thisgen=num2skim;
gp.state.newpopbuild_in_thisgen=num2build;
% disp(['pareto number in this generation:    ' num2str(gp.state.pareto_in_thisgen)]);
% disp(['elite number in this generation:    ' num2str(gp.state.elite_in_thisgen)]);
% disp(['new popbuild number in this generation:    ' num2str(gp.state.newpopbuild_in_thisgen)]);


%gp.runcontrol.pop_size = num2skim+num2build;


%popbuild_par(gp)



