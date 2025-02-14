function [check_index,newPop_out] = popbuild_eT3(tempgp)


max_nodes = tempgp.treedef.max_nodes;
maxNodesInf = isinf(max_nodes);
maxDepth = tempgp.treedef.max_depth;
useMultiGene = tempgp.genes.multigene;
p_cross_hi = tempgp.genes.operators.p_cross_hi;
crossRate = tempgp.genes.operators.hi_cross_rate;
max_number_of_ERC=tempgp.nodes.inputs.max_number_of_ERC;
maxGenes = tempgp.genes.max_genes;

highLevelCross = false;
if useMultiGene
    
    %select crossover type if multigene enabled
    if rand < p_cross_hi
        highLevelCross = true;
    end
    
end
%Select Parents
parentIndex = selection(tempgp);
dad = tempgp.pop{parentIndex};
numDadGenes = numel(dad);

parentIndex = selection(tempgp);
mum = tempgp.pop{parentIndex};
numMumGenes = numel(mum);

if highLevelCross
    if numMumGenes>1 || numDadGenes>1
        %high level crossover (only use if either parent has more
        %than one gene) This is modified/simplified from the
        %version in GPTIPS v1. This version just chooses between 1
        %and N genes to exchange randomly and independently for
        %each parent with a probability set by
        %gp.genes.operators.hi_cross_rate (where N =  number of
        %genes in that individual). Like GPTIPS 1, this assumes
        %that the order of genes within an individual is not
        %important.
        
        dadGeneSelectionInds = rand(1,numDadGenes) < crossRate;
        mumGeneSelectionInds = rand(1,numMumGenes) < crossRate;
        
        if ~any(dadGeneSelectionInds)
            dadGeneSelectionInds(1,ceil(numDadGenes *rand)) = true;
        end
        
        if ~any(mumGeneSelectionInds)
            mumGeneSelectionInds(1,ceil(numMumGenes *rand)) = true;
        end
        
        dadSelectedGenes = dad(dadGeneSelectionInds);
        mumSelectedGenes = mum(mumGeneSelectionInds);
        
        dadRemainingGenes = dad(~dadGeneSelectionInds);
        mumRemainingGenes = mum(~mumGeneSelectionInds);
        
        mumOffspring = [mumRemainingGenes dadSelectedGenes];
        dadOffspring = [dadRemainingGenes mumSelectedGenes];
        
        %curtail offspring longer than the max allowed number of genes
        %before writing to new population (only write 1 if no space for 2 offspring)
        %newPop{buildCount,1} = mumOffspring(1:(min(end,maxGenes)));
        %buildCount = buildCount+1;
        [newPop_temp1,check_index1] = orderconst(mumOffspring(1:(min(end,maxGenes))),max_number_of_ERC);
        [newPop_temp2,check_index2] = orderconst(dadOffspring(1:(min(end,maxGenes))),max_number_of_ERC);
        if (check_index1==1)&&(check_index2==1)
            newPop_out={newPop_temp1,newPop_temp2};
            check_index=1;
        elseif (check_index1==1)&&(check_index2==0)
            newPop_out={newPop_temp1,[]};
            check_index=1;
        elseif (check_index1==0)&&(check_index2==1)
            newPop_out={newPop_temp2,[]};
            check_index=1;
        else
            check_index=0;
            newPop_out=[];
        end
        
    else
        highLevelCross = false;
    end
end

%low level crossover: picks a random gene from each parent and
%crosses them. The offspring replace the original target genes
if ~highLevelCross
    
    if useMultiGene %if multigene then select target genes
        dad_target_gene_num = ceil(rand*numDadGenes); %randomly select a gene from dad
        mum_target_gene_num = ceil(rand*numMumGenes); %randomly select a gene from mum
        dad_target_gene = dad{1,dad_target_gene_num};
        mum_target_gene = mum{1,mum_target_gene_num};
    else
        dad_target_gene_num = 1;
        mum_target_gene_num = 1;
        dad_target_gene = dad{1};
        mum_target_gene = mum{1};
    end
    
    
    for loop = 1:10  %loop (max 10 times) until both children meet size constraints
        
        %produce 2 offspring
        [son,daughter] = crossover(mum_target_gene,dad_target_gene,tempgp);
        son_depth = getdepth(son);
        
        %check if both children meet size and depth constraints
        %if true then break out of while loop and proceed
        crossOverSuccess = false;
        if son_depth <= maxDepth
            daughter_depth = getdepth(daughter);
            if daughter_depth <= maxDepth
                
                if maxNodesInf
                    crossOverSuccess = true;
                    break;
                end
                
                son_nodes = getnumnodes(son);
                if son_nodes <= max_nodes
                    daughter_nodes = getnumnodes(daughter);
                    if  daughter_nodes <= max_nodes
                        crossOverSuccess = true;
                        break;
                    end
                end
            end
        end
        
    end
    
    %if no success then re-insert parents
    if ~crossOverSuccess
        son = dad_target_gene;
        daughter = mum_target_gene;
    end
    
    %write offspring back to right gene positions in parents and write to population
    dad{1,dad_target_gene_num} = son;
    mum{1,mum_target_gene_num} = daughter;
    [newPop_temp1,check_index1] = orderconst(dad,max_number_of_ERC);
    [newPop_temp2,check_index2] = orderconst(mum,max_number_of_ERC);
    if (check_index1==1)&&(check_index2==1)
        newPop_out={newPop_temp1,newPop_temp2};
        check_index=1;
    elseif (check_index1==1)&&(check_index2==0)
        newPop_out={newPop_temp1,[]};
        check_index=1;
    elseif (check_index1==0)&&(check_index2==1)
        newPop_out={newPop_temp2,[]};
        check_index=1;
    else
        check_index=0;
        newPop_out=[];
    end
   
end %end of if ~use_high


% savePath=[tempgp.fitness.savepath, '666.mat'];
% save(savePath);








end

