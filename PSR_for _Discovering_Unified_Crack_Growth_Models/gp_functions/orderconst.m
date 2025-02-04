function [newPop_out,check_index] = orderconst(newPop_in,max_const_number)
%ORDERCONST order const c of new pop
     newPop_temp=newPop_in;
     numGenes = numel(newPop_temp);
     Const_new_number=0;
     for i=1:numGenes
         posC=strfind(newPop_temp{i},'c');
         for j=1:numel(posC)
            Const_new_number=Const_new_number+1;
            newPop_temp{i} = strrep(newPop_temp{i},newPop_temp{i}(posC(j)+1),sprintf('%d',Const_new_number));
         end
     end
     if Const_new_number<max_const_number+1
        newPop_out=newPop_temp;
        check_index=1;
     else
        newPop_out=newPop_temp;
        check_index=0;
     end
end

