function [ypred,theta] = ypredtrainvalues( gp,i,j )

if nargin<3
    j=1;
end
    evalstr = tree2evalstr(gp.pop{i},gp);
    xtrain=[gp.multidata{1,j}(:,1),gp.multidata{1,j}(:,4)];
    ytrain=log10(gp.multidata{1,j}(:,2));
    % process evalstr with regex to allow direct access to data matrices
    pat = 'x(\d+)';
    evalstr = regexprep(evalstr,pat,'xtrain(:,$1)');
    y = ytrain;
    [numData,~] =size(ytrain);
    numGenes = numel(evalstr); 
    %set up a matrix to store the tree outputs plus a bias column of ones
    geneOutputs = ones(numData,numGenes+1);
    
    %eval each gene in the current individual
    for i = 1:numGenes
        ind = i + 1;
        eval(['geneOutputs(:,ind)=' evalstr{i} ';']);
    end
    goptrans = geneOutputs';
    prj = goptrans * geneOutputs;
    theta = pinv(prj) * goptrans * y;
    ypred = geneOutputs * theta;

    


end


