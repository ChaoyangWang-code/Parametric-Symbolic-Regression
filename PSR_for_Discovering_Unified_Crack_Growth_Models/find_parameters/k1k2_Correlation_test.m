function [k1k2_Correlation] = k1k2_Correlation_test(evalstr_IN)
%K1K2_CORRELATION_TEST 
numGenes = numel(evalstr_IN); 
num_const=0;
for i=1:numGenes
    open_sq_br = strfind(evalstr_IN{i},'c');
    num_const = num_const+numel(open_sq_br);
end

for i=1:numGenes+2
    syms(['f',num2str(i)]);
end
for i=1:num_const
    syms(['c',num2str(i)]);
end


evalstr1=evalstr_IN;
evalstr2=evalstr_IN;

syms delta_K Kmax k1 k2
%eq=['log10(f1)'];

x1_str=['(','delta_K/','k1',')'];
x2_str=['(','Kmax/','k2',')'];    
for i=1:numGenes
    evalstr1{i} = strrep(evalstr1{i},'x1',x1_str);
    evalstr1{i} = strrep(evalstr1{i},'x2',x2_str);
    %eq_part{i}=['f',num2str(i+1),'*',evalstr_IN{i}];
    eq_part1{i}=['f',num2str(i+1),'*','log10(',evalstr1{i},')',];
	if i==1
		eq1=eq_part1{i};
	else
		eq1=[eq1,'+',eq_part1{i}];
	end
end


x1_str=['(','delta_K/','(pi*k1)',')'];
x2_str=['(','Kmax/','(pi*k2)',')'];    
for i=1:numGenes
    evalstr2{i} = strrep(evalstr2{i},'x1',x1_str);
    evalstr2{i} = strrep(evalstr2{i},'x2',x2_str);
    %eq_part{i}=['f',num2str(i+1),'*',evalstr_IN{i}];
    eq_part2{i}=['f',num2str(i+1),'*','log10(',evalstr2{i},')',];   
	if i==1
		eq2=eq_part2{i};
	else
		eq2=[eq2,'+',eq_part2{i}];
	end	
end


warning('off', 'all');

try
    expr1=simplify(eval(eq1));
    expr2=simplify(eval(eq2));
    if isequal(expr1,expr2)
        k1k2_Correlation=1;
    elseif isAlways(eval(eq1) == eval(eq2))
        k1k2_Correlation=1;
    else
        k1k2_Correlation=0;
    end
catch
    try
        if isAlways(eval(eq1) == eval(eq2))
            k1k2_Correlation=1;
        else
            k1k2_Correlation=0;
        end
    catch
            k1k2_Correlation=0;
    end

end

%k1k2_Correlation=isequal(simplify(eval(eq1)),simplify(eval(eq2)));
%k1k2_Correlation=isAlways(eval(eq1) == eval(eq2));


warning('on', 'all');
end
%k1k2_Correlation_test(evalstr_in)
