function [ out ] = constrain( in, lower, upper )
% constrains the input between the upper and lower bounds

in(in < lower) = lower;
in(in > upper) = upper;

in(isnan(in)) = lower;

out = in;


end

