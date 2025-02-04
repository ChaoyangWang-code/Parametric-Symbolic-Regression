function [new_pop_d] = delete_same_gene(old_pop,num_genes,tempgp)
    syms x1 x2 c1 c2       

    old_pop_evalstr=tree2evalstr(old_pop,tempgp);             
    new_pop_index=1;
    %new_pop_d = cell(1,1); % 
    %new_pop_d_evalstr = cell(1,1);
    warning('off', 'all');
    for i = 1:num_genes
        if i==1
            new_pop_d{1,1}=old_pop{1,1};
            new_pop_d_evalstr{1,1}=old_pop_evalstr{1,1};
            new_pop_index=new_pop_index+1;
        else
            is_duplicate = false; % 
            for j = 1:size(new_pop_d,2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                try
                    expr1=simplify(eval(old_pop_evalstr{1,i}));
                    expr2=simplify(eval(new_pop_d_evalstr{1,j}));
                    if isequal(expr1,expr2)
                        is_duplicate = true; 
                        break; 
                    elseif isAlways(eval(old_pop_evalstr{1,i}) == eval(new_pop_d_evalstr{1,j}))
                        is_duplicate = true; 
                        break; 
                    end
                catch
                    try
                        if isAlways(eval(old_pop_evalstr{1,i}) == eval(new_pop_d_evalstr{1,j}))
                            is_duplicate = true; 
                            break; 
                        end
                    catch
                            %k1k2_Correlation=0;
                    end
                
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if isequal(simplify(eval(old_pop_evalstr{1,i})),simplify(eval(new_pop_d_evalstr{1,j})))      
                %     is_duplicate = true; 
                %     break; 
                % elseif isAlways(eval(old_pop_evalstr{1,i}) == eval(new_pop_d_evalstr{1,j}))
                %     is_duplicate = true; 
                %     break; 
                % end
            end
            if ~is_duplicate
                new_pop_d{1,new_pop_index} = old_pop{1,i}; 
                new_pop_d_evalstr{1,new_pop_index}=old_pop_evalstr{1,i};
                new_pop_index=new_pop_index+1;
            end
        end
    end
    warning('on', 'all');
%delete_same_gene(newPop{i,1},3)


end

