
clc
clear all
close all

color = [
    1, 0, 0;         
    0, 1, 0;         
    0, 0, 1;         
    1, 0.843, 0;     
    0.580, 0, 0.827; 
    1, 0.647, 0;     
    0.529, 0.808, 0.922; 
    0.133, 0.545, 0.133; 
    1, 0.078, 0.576; 
    0, 0.808, 0.820  
];

%figuresize=[100, 100, 1200, 1000]; %[left bottom width height]
for gen_count=60
    figuresize=[100, 100, 1000, 800]; %[left bottom width height]
    %color = [1 0.38 0;1 0 0; 0 0 1;0.63 0.13 0.94;0.5 1 0;1 0.5 0.31; 0.53 0.15 0.34; 1 0 1;1 1 0; 0.45 0.29 0.07;0 0 1;0.03 0.18 0.33];  
    colors=[color;color;color;color];
    Path='workdir';     %Select the result path,example:    Windows:'C:\Users\Desktop\test250106\'   linux:'/gs/home/gp_test/test1/';
    show_figure='on';  
    str=num2str(gen_count);
    data_path=[Path, str,'\gen',str, '.mat'];   
    data_geni{gen_count}=load(data_path);
    
    
    foldername_gen=[Path, str];
    %show_results_gen_i(figuresize,gp,foldername_gen,str)
    best_pop=MULTI_R_show_results_of_all_pareto_in_gen_i(colors,show_figure,figuresize,data_geni{gen_count}.gp,foldername_gen,str);
    clc
end


function [ best_pop ] = MULTI_R_show_results_of_all_pareto_in_gen_i(colors,show_figure,figuresize,gp,foldername_gen,str)
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
        error('All inf');
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
    %% pareto
    if gp.fitness.gen_count_now<=gp.runcontrol.stage1
        h0 = figure('Visible', show_figure);
        hold on;grid on
        title('pareto pops')
        plot(pareto_comp,pareto_fit,'o', 'MarkerFaceColor', 'b')        
    elseif (gp.fitness.gen_count_now>gp.runcontrol.stage1)&&(gp.fitness.gen_count_now<=gp.runcontrol.stage2)
        h0 = figure('Visible', show_figure);
        hold on;grid on
        title('pareto pops')
        plot(pareto_comp,pareto_fit,'o', 'MarkerFaceColor', 'b')              
    else
        h0 = figure('Visible', show_figure);
        hold on;grid on
        title('pareto pops')
        plot(pareto_comp,pareto_fit,'o', 'MarkerFaceColor', 'b')                
    end

    %% 
    returns=gp.fitness.returnvalues;
    for pareto_i=1:pu_num    
        str_pareto=num2str(pareto_i);
        pvals=rank1_tour_ind(pareto_i);
        pop_now=returns{pvals};
        parameter_now=pop_now{1};
        numGenes=pop_now{2};
        num_const=pop_now{3};    
        F_eq=pop_now{5};
        for i=1:numGenes+2
            syms(['f',num2str(i)]);
        end
        for i=1:num_const
            syms(['c',num2str(i)]);
        end
        
        F_eq_show_temp=F_eq;
        syms delta_K Kmax k1 k2
        eq=['log10(f1)'];

        for i=1:numGenes
            x1_str=['(','delta_K/','k1',')'];
            x2_str=['(','Kmax/','k2',')'];
            F_eq_show_temp{i} = strrep(F_eq_show_temp{i},'x1',x1_str);
            F_eq_show_temp{i} = strrep(F_eq_show_temp{i},'x2',x2_str);
            %eq_part{i}=['f',num2str(i+1),'*',F_eq_show_temp{i}];
            eq_part{i}=['f',num2str(i+1),'*','log10(',F_eq_show_temp{i},')',];
            eq=[eq,'+',eq_part{i}];
        end
        eq=[eq,'+','f',num2str(numGenes+2),'*','log10(delta_K)'];
        eq10=['10^(',eq,')'];

        pat2 = 'c(\d+)';
        Const_pair_now=parameter_now(1,4:3+num_const);   

        F_eq_show=eq;
        F_eq_show10=eq10;

        materials_N=size(gp.multidata,2);                  
        data_N=size(gp.all_data_without_collection,2);     
        data_all=gp.all_data_without_collection;
        data_all_collection=gp.multidata;
        foldername_pareto_i_fig=[foldername_gen,'\','pareto_',num2str(pareto_i),'_fig'];
        if ~exist(foldername_pareto_i_fig, 'dir')
            mkdir(foldername_pareto_i_fig);
        end        
        foldername_pareto_i_png=[foldername_gen,'\','pareto_',num2str(pareto_i),'_png'];
        if ~exist(foldername_pareto_i_png, 'dir')
            mkdir(foldername_pareto_i_png);
        end  
        %% 
        for i=1:materials_N
            
            parameter_K=parameter_now(i,2:3);
            delta_kth=parameter_K(1);
            kc=parameter_K(2);
            Const_pair_now=parameter_now(i,4:3+num_const);
            theta=parameter_now(i,4+num_const:4+num_const+numGenes+1);           
            %
            D_i=gp.multidata{1,i};
            %          
            material_i_delta_K_min=min(data_all_collection{1,i}(:,1));
            %material_i_delta_K_min=6;
            material_i_delta_K_max=max(data_all_collection{1,i}(:,1));
            material_i_R_min=min(data_all_collection{1,i}(:,3));
            material_i_R_max=max(data_all_collection{1,i}(:,3));            
            % material_i_delta_K_pre=0.9*material_i_delta_K_min:0.01:1.1*material_i_delta_K_max;
            % material_i_delta_K_pre=material_i_delta_K_pre';
            material_i_delta_K_pre = linspace(material_i_delta_K_min, material_i_delta_K_max,100);
            material_i_R_pre = linspace(material_i_R_min, material_i_R_max,100);  
            [material_i_delta_K_PREDICTION, material_i_R_PREDICTION] = meshgrid(material_i_delta_K_pre, material_i_R_pre);
            xtest{1}=material_i_delta_K_PREDICTION./delta_kth;
            xtest{2}= material_i_delta_K_PREDICTION./(1-material_i_R_PREDICTION)./kc;
            %xtrain=[material_i_delta_K_P./delta_kth,   material_i_delta_K_P./(1-material_i_R_P)./kc];
            num_of_predata_i=size(material_i_delta_K_pre,2);
            pat1 = 'x(\d+)';
            pat2 = 'c(\d+)';
            F_eq=pop_now{5};
            F_eq = regexprep(F_eq,pat2,'Const_pair_now($1)');
            F_eq = regexprep(F_eq,pat1,'xtest{$1}');                                         
            geneOutputs_test = theta(1).*ones(num_of_predata_i,num_of_predata_i);
            for j = 1:numGenes
                gene_temp=eval([F_eq{j} ';']);
                geneOutputs_test=geneOutputs_test+theta(j+1).*lg(gene_temp);
            end 
            geneOutputs_test=geneOutputs_test+theta(2+numGenes).*lg(material_i_delta_K_PREDICTION);       
            %nothasComplex = isreal(geneOutputs); 
            if isreal(geneOutputs_test) && ~any(isinf(geneOutputs_test),'all')
                test_index = 1;
            else
                test_index = 0; 
            end

            %if test_index
            try 
                figure('Visible', show_figure)            
                %%%%%%%
                scatter3(log10(D_i(:,1)),D_i(:,3),log10(D_i(:,2)),'o','MarkerEdgeColor', colors(i, :),'MarkerFaceColor', colors(i, :))
                hold on;
                %%%%%%
                h1 = surf(lg(material_i_delta_K_PREDICTION), material_i_R_PREDICTION, geneOutputs_test);
              
                set(h1, 'EdgeColor', 'none');            
                set(h1, 'FaceColor', 'blue');
                set(h1, 'FaceAlpha', 0.3);
                set(gca, 'YDir', 'reverse')
                %xlabel('$\Delta \mathit{K}\ (\mathrm{ksi\sqrt{in}})$', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
                %ylabel('$\mathit{R}$', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
                %zlabel('$\mathrm{d}\mathit{a}/\mathrm{d}\mathit{N}\ (\mathrm{in/cycle})$', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
                %title(['Materials-', num2str(i)],'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 14,'FontWeight', 'bold');
                % title(['$\textit{Materials-', num2str(i), '}$'], 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
                %title(['Data label:',num2str(i)], 'FontName', 'Arial', 'FontSize', 20, 'FontWeight', 'bold')
                %title(['pop:',num2str(pvals), ':' , '$\textit{Materials-', num2str(i), '}$'],'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold')
                set(gca, 'FontName', 'Arial', 'FontSize', 24);

                if min(geneOutputs_test(:))<-8
                    zlim([min(geneOutputs_test(:)) max(geneOutputs_test(:))]); 
                else
                    zlim([log10(1e-8) max(geneOutputs_test(:))]); 
                end

                %zlim([log10(1e-10) log10(1)]); 
                % z_ticks = log10([1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1]);
                % z_labels = {'10^{-10}', '10^{-9}', '10^{-8}','10^{-7}', '10^{-6}', '10^{-5}','10^{-4}', '10^{-3}', '10^{-2}','10^{-1}'};
                z_ticks = log10([1e-10 1e-8 1e-6 1e-4 1e-2 ]);
                z_labels = {'10^{-10}', '10^{-8}', '10^{-6}', '10^{-4}', '10^{-2}'};                
                set(gca, 'ZTick', z_ticks, 'ZTickLabel', z_labels);
                if material_i_delta_K_max>100
                    xlim([log10(1) log10(material_i_delta_K_max)]); 
                else
                    xlim([log10(1) log10(100)]); 
                end
                
                % x_ticks = log10([1 10 20 30 40 80 120]);
                % x_labels = {'1', '10', '20','30','40','80','120'};
                x_ticks = log10([1 10 100]);
                x_labels = {'1', '10', '100'};               
                set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels);

            %else
            catch
                xtrain=[D_i(:,1)./delta_kth,D_i(:,4)./kc];
                num_of_data_i=size(D_i,1);
                F_eq=pop_now{5};
                pat1 = 'x(\d+)';
                pat2 = 'c(\d+)';
                F_eq = regexprep(F_eq,pat2,'Const_pair_now($1)');
                F_eq = regexprep(F_eq,pat1,'xtrain(:,$1)');
                geneOutputs = ones(num_of_data_i,numGenes+2);
                for j = 1:numGenes
                    ind = j + 1;
                    gene_temp=eval_exe(F_eq{j},xtrain,Const_pair_now);
                    geneOutputs(:,ind)=lg(gene_temp);
                end
                geneOutputs(:,numGenes+2)=lg(D_i(:,1));                 
                ypredtrain = geneOutputs * theta';

                h1 = figure('Visible', show_figure); 
                plot(lg(D_i(:,1)), lg(D_i(:,2)),'+','Color', colors(1, :),'MarkerSize', 10);
                hold on; grid on
                plot(lg(D_i(:,1)),ypredtrain,'o','Color', colors(3, :),'LineWidth',1.5)  
                %xlabel('$\Delta K\ (MPa\sqrt{m})$', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 14,'FontWeight', 'bold');
                %ylabel('$da/dN\ (mm/cycle)$', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 14,'FontWeight', 'bold');
                set(gca, 'FontName', 'Arial', 'FontSize', 12);     
            end
        end
        %% 
        theta_temp=parameter_now(:,4+num_const:4+num_const+numGenes+1); 
        theta = [theta_temp(:, 1), theta_temp(:, end), theta_temp(:, 2:end-1)];
        % 
        num_coefficients = size(theta, 2);
        % 
        coeff_labels = arrayfun(@(x) ['f' num2str(x)], 1:num_coefficients, 'UniformOutput', false);
        h2 = figure('Visible', show_figure);
        %BOX = boxplot(theta, 'Labels', coeff_labels);
        BOX = boxplot(theta);
        set(BOX, 'LineWidth', 2); 

        set(gca, 'FontName', 'Arial', 'FontSize', 20);
        %xlabel('Parameters', 'FontName', 'Arial', 'FontSize', 20);
        %ylabel('Value', 'FontName', 'Arial', 'FontSize', 20);
        %title(['Distribution of Coefficients for all labeled data'], 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
  
        xticks(1:num_coefficients);
        xticklabels(coeff_labels);
        xticklabels([]);

        %set(gca, 'Box', 'on', 'LineWidth', 2);
        grid on;
        set(gca, 'GridLineStyle', '-'); 
        set(gca, 'GridColor', [0.6 0.6 0.6]); 
        set(gca, 'FontSize', 20);

    end      

end

end














