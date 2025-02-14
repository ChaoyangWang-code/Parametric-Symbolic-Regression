function [check_index,newPop_out] = popbuild_eT1(parent,tempgp)

maxDepth = tempgp.treedef.max_depth;
max_nodes = tempgp.treedef.max_nodes;
maxNodesInf = isinf(max_nodes);

    if tempgp.genes.multigene %if using multigene, extract a target gene
        numParentGenes = numel(parent);
        targetGeneIndex = ceil(rand * numParentGenes); %randomly select a gene from parent
        targetGene = parent{1,targetGeneIndex}; %extract it
    else
        targetGeneIndex = 1;
        targetGene = parent{1};
    end

    mutateSuccess = false;
    for loop = 1:10	%loop until a successful mutation occurs (max loops=10)
        
        mutatedGene = mutate(targetGene,tempgp);
        mutatedGeneDepth = getdepth(mutatedGene);
        
        if mutatedGeneDepth <= maxDepth
            
            if maxNodesInf
                mutateSuccess = true;
                break;
            end
            
            mutatedGeneNodes = getnumnodes(mutatedGene);
            if mutatedGeneNodes <= max_nodes
                mutateSuccess = true;
                break;
            end
        end %end of constraint check
        
    end  %end of mutate for loop

    %if no success then use parent gene
    if ~mutateSuccess
        mutatedGene = targetGene;
    end

    %add the mutated individual to new pop
    parent{1,targetGeneIndex} = mutatedGene;
    %newPop{buildCount,1} = parent;
    [newPop_temp,check_index] = orderconst(parent,tempgp.nodes.inputs.max_number_of_ERC);
    if check_index==1
        newPop_out=newPop_temp;
    else
        newPop_out=[];
    end


end

