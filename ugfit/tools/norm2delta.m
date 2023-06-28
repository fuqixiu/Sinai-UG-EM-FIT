function [ out ] = norm2delta(in)
% transformation from gaussian space to alpha space, as suggested in Daw 2009 Tutorial
% BRK Shevlin 2023

out = 4*(logsig(in))-2; 

end

