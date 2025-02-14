function b = gpnot(a)
%GPNOT Node. Wrapper for logical NOT 
%
%   Copyright (c) 2009-2015 Dominic Searson 
%   Copyright (c) 2023-2025 Chaoyang Wang 
%   GPTIPS 2
%
%   See also GPOR, GPAND, LTH, GTH, MAXX, MINX, NEG, IFLTE, STEP, THRESH

if any(isnan(a))
    b = nan;
    return;
end

b = double(not(real(a)));