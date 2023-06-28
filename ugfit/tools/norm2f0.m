function [ out ] = norm2f0(in)
% transformation from gaussian space to alpha space, as suggested in Daw 2009 Tutorial
% BRK Shevlin 2023

out = 20*(logsig(in)); 

end

