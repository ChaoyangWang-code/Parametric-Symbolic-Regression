function gp=initbuild_copy(gp)
%INITBUILD Generate an initial population of GP individuals.
%
%   GP = INITBUILD(GP) creates an initial population using the parameters
%   in the structure GP. Various changes to fields of GP are made.
%
%   Remarks:
%   
%   Each individual in the population can contain 1 or more separate trees
%   (if more than 1 then each tree is referred to as a gene of the
%   individual). This is set by the user in the config file in the field
%   GP.GENES.MULTIGENE.
%
%   Each individual is a cell array. You can use cell addressing to
%   retrieve the individual genes. E.g. GP.POP{3}{2} will return the 2nd
%   gene in the third individual.
%
%   Trees are (by default) constructed using a probabilistic version of
%   Koza's ramped 1/2 and 1/2 method. E.g. if maximum tree depth is 5 and
%   population is 100 then, on average, 20 will be generated at each depth
%   (1/2 using 'grow' and 1/2 using 'full').
%
%   Multiple copies of genes in an individual are disallowed when
%   individuals are created. There is, however, no such restriction imposed
%   on future generations.
%
%   (c) Dominic Searson 2009-2015
%   Copyright (c) 2023-2025 Chaoyang Wang
%   GPTIPS 2
%
%   See also POPBUILD, TREEGEN

%load('pop_231116_1.mat');
%load('pop_test231129.mat');
%load('pop_231206.mat');                       
load('pop_test231214.mat');
popSize = gp.runcontrol.pop_size;



gp.pop = cell(popSize,1);
gp.pop =pop_0(1:popSize,1);
clear pop_0;

