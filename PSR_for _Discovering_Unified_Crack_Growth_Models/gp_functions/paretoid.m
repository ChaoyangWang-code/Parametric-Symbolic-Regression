function paretoID = paretoid(gp)

geshu = 1;
numpop=gp.runcontrol.pop_size;
paretoID=[];
fitness_all(:,1)=gp.fitness.values;
fitness_all(:,2)=gp.fitness.complexity;
for i=1:numpop
    pareto_selection=true;
    beizhipeiyinzi=false;
    if fitness_all(i,1)==inf
        pareto_selection=false;
    else
        if i==1
            for j=2:numpop
                if iszhipei(i,j,fitness_all)
                pareto_selection=false;
                 break
                end
            end 
        elseif i~=numpop
            for j=1:(i-1)
                 if iszhipei(i,j,fitness_all)  || determinate_equal(gp.pop{i,1},gp.pop{j,1})
                     pareto_selection=false;
                     beizhipeiyinzi=true;
                     break
                 end
            end
            if beizhipeiyinzi==false
             for j=(i+1):numpop
                 if iszhipei(i,j,fitness_all)
                     pareto_selection=false;
                     break
                 end
             end
            end
        else
            for j=1:numpop-1
                if iszhipei(i,j,fitness_all)||  determinate_equal(gp.pop{i,1},gp.pop{j,1})
                pareto_selection=false;
                break
                end
            end
        end
    end
    if pareto_selection==true
        paretoID(geshu,1)=i;
        geshu=geshu+1;
    end
end
fitnessall_pareto=fitness_all(paretoID,:);
paretoID=[paretoID,fitnessall_pareto];
end

function determination =determinate_equal(a,b)
%
e=[];
determination = false;
if numel(a)~= numel(b)
    determination=false;
    return
else
    num = numel(a);
    num_equal = 0;
    for c = 1:num
        f = 0;
        for d = 1:num
            if isequal(a{1,c},b{1,d})
                if isempty(find(e==d))
                    num_equal = num_equal+1;
                    e = [e,d];
                    f=1;
                    break
                end
            end
        end
        if f == 0
            determination = false;
            return
        end
    end
    if num_equal == num
        determination = true;
    end
end
end





%
function a = iszhipei(i,j,fitness_all)
b=0;
if fitness_all(i,1)>fitness_all(j,1)
    b=b-1;
elseif fitness_all(i,1)<fitness_all(j,1)
    b=b+1;
end
if fitness_all(i,2)>fitness_all(j,2)
    b=b-1;
elseif fitness_all(i,2)<fitness_all(j,2)
    b=b+1;
end
if b<0
    a=1;
else
    a=0;
end
end
    
        
            
    
    