function [pareto_this_gen] = Find_pareto_pop_in_gen_i(gp)

if gp.fitness.loss_mode>=2
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
    num_rank1 = numel(rank1_tour_ind_temp);
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
    for pareto_i=1:pu_num     
        str_pareto=num2str(pareto_i);
        pvals=rank1_tour_ind(pareto_i);
        pareto_this_gen{pareto_i,1}=pvals;
        pareto_this_gen{pareto_i,2}=gp.fitness.returnvalues{pvals};
        pareto_this_gen{pareto_i,3}=gp.fitness.values_alldata{pvals};
        pareto_this_gen{pareto_i,4}=gp.fitness.values(pvals);
        pareto_this_gen{pareto_i,5}=gp.fitness.complexity(pvals);

    end

end


end

