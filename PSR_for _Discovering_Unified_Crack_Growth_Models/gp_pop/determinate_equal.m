function determination =determinate_equal(a,b)

e=[];   
determination = false;
if numel(a) ~= numel(b)
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
                    f = 1;
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