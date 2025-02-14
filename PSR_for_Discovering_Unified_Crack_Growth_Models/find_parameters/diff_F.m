function [eq,diff_omega,contains_k1,contains_k2 ] = diff_F(evalstr_IN,numConst,diff_model)

numGenes = numel(evalstr_IN); 
for i=1:numGenes+2
    syms(['f',num2str(i)]);
end
for i=1:numConst
    syms(['c',num2str(i)]);
end

if diff_model==2.1
    syms delta_K Kmax k1 k2
    %eq=['log10(f1)'];
    eq=['f1'];
    for i=1:numGenes
        x1_str=['(','delta_K/','k1',')'];
        x2_str=['(','Kmax/','k2',')'];
        evalstr_IN{i} = strrep(evalstr_IN{i},'x1',x1_str);
        evalstr_IN{i} = strrep(evalstr_IN{i},'x2',x2_str);
        %eq_part{i}=['f',num2str(i+1),'*',evalstr_IN{i}];
        eq_part{i}=['f',num2str(i+1),'*','log10(',evalstr_IN{i},')',];
        eq=[eq,'+',eq_part{i}];
    end
    eq=[eq,'+','f',num2str(numGenes+2),'*','log10(delta_K)'];
    variables_in_eq = symvar(eval(eq));
    % Check if k1 k2 exists.
    contains_k1 = ismember('k1', variables_in_eq);
    contains_k2 = ismember('k2', variables_in_eq);
%
    for i=1:numGenes+2
        diff_omega{i}=char(diff(eval(eq),['f',num2str(i)]));     
    end

    diff_omega{numGenes+3}=char(diff(eval(eq),k1));
    diff_omega{numGenes+4}=char(diff(eval(eq),k2));
    

    diff_omega=strrep(diff_omega,'/','./');
    diff_omega=strrep(diff_omega,'^','.^');
    diff_omega=strrep(diff_omega,'*','.*');

    eq=strrep(eq,'/','./');
    eq=strrep(eq,'^','.^');
    eq=strrep(eq,'*','.*');


elseif diff_model==2.2
    syms z1 z2 R pk1 pk2
    eq=['f',num2str(numGenes+2),'*','z1'];
    for i=1:numGenes
        x1_str=['(','pk1','/(10^z1)',')'];
        x2_str=['(','(10^z1)/','pk1',')'];
        x3_str=['(','(10^z2)/','pk2',')'];
        x4_str=['(','pk2','/(10^z2)',')'];
        x5_str='R';
        evalstr_IN{i} = strrep(evalstr_IN{i},'x1',x1_str);
        evalstr_IN{i} = strrep(evalstr_IN{i},'x2',x2_str);
        evalstr_IN{i} = strrep(evalstr_IN{i},'x3',x3_str);
        evalstr_IN{i} = strrep(evalstr_IN{i},'x4',x4_str);
        evalstr_IN{i} = strrep(evalstr_IN{i},'x5',x5_str);
        eq_part{i}=['f',num2str(i+1),'*',evalstr_IN{i}];
        eq=[eq,'+',eq_part{i}];
    end
    eq=[eq,'+','f1'];
    dF_dk1=char(diff(eval(eq),z1));                %dy/dk1
    dF_dk2=char(diff(eval(eq),z2));                %dy/dk2
    

    dF_dk1=strrep(dF_dk1,'/','./');
    dF_dk2=strrep(dF_dk2,'/','./');
    dF_dk1=strrep(dF_dk1,'^','.^');
    dF_dk2=strrep(dF_dk2,'^','.^');
    dF_dk1=strrep(dF_dk1,'*','.*');
    dF_dk2=strrep(dF_dk2,'*','.*');    

end

end

