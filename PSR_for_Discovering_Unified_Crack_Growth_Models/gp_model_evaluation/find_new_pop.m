function gp= find_new_pop(gp)
%FIND_NEW_POP 
    
    pop_in_this_gen=gp.pop;
    pop_size_this_gen=gp.runcontrol.pop_size;
    pop_all_history=gp.result_all_history(:,1);
    pop_size_all_history=size(pop_all_history,1);
    positions  = zeros( pop_size_this_gen,1); 


    parfor i = 1:pop_size_this_gen
        for j = 1:pop_size_all_history
            if determinate_equal(pop_in_this_gen{i}, pop_all_history{j}) 
                positions (i,1) = j;       
                break; % 
            end
        end
    end

    toc

    matchedIndices = positions(positions > 0);  
    B1 = pop_all_history(matchedIndices);     
    B2 = pop_in_this_gen(positions == 0);  
    gp.state.new_pop_num=size(B2,1);
    % pop_size_this_gen_new=[B2,B1];
    % gp.pop=pop_size_this_gen_new;
    pop_size_this_gen_new=B2;
    gp.pop=pop_size_this_gen_new;    
    gp.positions=matchedIndices;
end

%find_new_pop(pop_all,gp)




