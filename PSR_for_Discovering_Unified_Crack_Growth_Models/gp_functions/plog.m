function out = plog(x)
%PLOG Node. Calculate the element by element protected natural log of a vector. 
%
%   OUT = PLOG(X) performs LOG(ABS(X))
%
%   Copyright (c) 2009-2015 Dominic Searson
%   Copyright (c) 2023-2025 Chaoyang Wang
%   GPTIPS 2
%
%   See also PDIV

out = log(abs(x));
out(isinf(out)) = 0;
