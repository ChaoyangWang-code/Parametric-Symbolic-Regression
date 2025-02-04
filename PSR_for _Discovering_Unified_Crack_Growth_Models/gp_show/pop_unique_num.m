function [uniqueCount ] = pop_unique_num(gp)
% Check how many unique pops there are.
    pop_all=gp.pop;
    pop_size=gp.runcontrol.pop_size;
    isUnique = true(1, pop_size); 


    parfor i = 1:pop_size
        for j = i+1:pop_size
            if isUnique(i) && (determinate_equal(pop_all{i}, pop_all{j}) )%
                isUnique(i) = false; 
                break; 
            end
        end
    end

uniqueCount  = sum(isUnique);


end

