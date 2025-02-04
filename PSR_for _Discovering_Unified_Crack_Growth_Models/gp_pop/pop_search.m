
for gen_count=1:16


    Path='B:\Desktop\gp_plot\';
    show_figure='on';
    
    str=num2str(gen_count);
    
    data_path=[Path, str,'\gen',str, '.mat'];
    
    data_geni{gen_count}=load(data_path);
    
    gp=data_geni{gen_count}.gp;
    
    returns=gp.fitness.returnvalues;
    
    for i=1:gp.runcontrol.pop_size
        pop_now=returns{i};        
        numGenes=pop_now{2};
        indi=pop_now{1,4};
        if numGenes==2
            if strcmp(indi{1,1},'minus(c1,f_r(x1))')&&strcmp(indi{1,2},'minus(c2,x2)')
                disp('bingo');
            elseif strcmp(indi{1,1},'minus(c1,x2)')&&strcmp(indi{1,2},'minus(c2,f_r(x1))')
                disp('bingo');
            end
        end
    
    end


end






