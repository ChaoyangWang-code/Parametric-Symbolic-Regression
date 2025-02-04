function y = psqroot(x)
%PSQROOT Node. Computes element by element protected square root of a vector
%
%   Y = PSQROOT(X) performs the vector operation Y = SQRT(ABS(X))   
%
%   Copyright (c) 2009-2015 Dominic Searson 
%   Copyright (c) 2023-2025 Chaoyang Wang  
%   GPTIPS 2
%
%   See also SQUARE, CUBE
y = sqrt(abs(x));
