function [ qout,qbounds ] = norm2par(modelID,qin)
% Lookup table to transform gaussian parameters to real parameter values for all models;
%  
% columns are different free parameters, have to be in the order defined by the model
% Created by MKW 2017
% Adapted by BRK Shevlin 2023
%
% INPUT:    - qin: input parameters, observation x parameter
% OUTPUT:   - qout: qin, transformed
%           - qbounds: sensible bounds for the parameter, rows are min/max, cols are parameters

%% in case input is a vector

if size(qin,2)==1, qin = qin'; end                                           % transpose input if it has wrong dimension, but only if it is a n x 1 vector


%% main models
% norm2alpha uses sigmoid (0 to 1)
% norm2beta uses exponential
% norm2sblur uses 

if strcmp(modelID, 'ms_UG0_adaptiveNorm')
    if size(qin,2)~=4, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2f0(qin(:,3)) norm2alpha(qin(:,4))];
    qbounds = [0 1; 0 15; 0 20; 0 1]';
elseif strcmp(modelID, 'ms_UG0_f0f_adaptiveNorm')
    if size(qin,2)~=3, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2alpha(qin(:,3))];
    qbounds = [0 1; 0 15; 0 1]';
elseif strcmp(modelID, 'ms_UG1_etaf_adaptiveNorm')
    if size(qin,2)~=5, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2f0(qin(:,3)) norm2alpha(qin(:,4)) norm2delta(qin(:,5))];
    qbounds = [0 1; 0 15; 0 20; 0 1; -2 2]';
elseif strcmp(modelID, 'ms_UG1_etaf_f0f_adaptiveNorm')
    if size(qin,2)~=4, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2alpha(qin(:,3)) norm2delta(qin(:,4))];
    qbounds = [0 1; 0 15; 0 1; -2 2]';
elseif strcmp(modelID, 'ms_UG2_etaf_adaptiveNorm')
    if size(qin,2)~=5, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2f0(qin(:,3)) norm2alpha(qin(:,4)) norm2delta(qin(:,5))];
    qbounds = [0 1; 0 15; 0 20; 0 1; -2 2]';
elseif strcmp(modelID, 'ms_UG2_etaf_f0f_adaptiveNorm')
    if size(qin,2)~=4, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2alpha(qin(:,3)) norm2delta(qin(:,4))];
    qbounds = [0 1; 0 15; 0 1; -2 2]';
elseif strcmp(modelID, 'ms_UG3_etaf_adaptiveNorm')
    if size(qin,2)~=5, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2f0(qin(:,3)) norm2alpha(qin(:,4)) norm2delta(qin(:,5))];
    qbounds = [0 1; 0 15; 0 20; 0 1; -2 2]';
elseif strcmp(modelID, 'ms_UG3_etaf_f0f_adaptiveNorm')
    if size(qin,2)~=4, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2alpha(qin(:,3)) norm2delta(qin(:,4))];
    qbounds = [0 1; 0 15; 0 1; -2 2]';
elseif strcmp(modelID, 'ms_UG0_fixedNorm')
    if size(qin,2)~=3, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2f0(qin(:,3))];
    qbounds = [0 1; 0 15; 0 20]';
elseif strcmp(modelID, 'ms_UG1_etaf_fixedNorm')
    if size(qin,2)~=4, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2f0(qin(:,3)) norm2delta(qin(:,4))];
    qbounds = [0 1; 0 15; 0 20; -2 2]';
elseif strcmp(modelID, 'ms_UG2_etaf_fixedNorm')
    if size(qin,2)~=4, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2f0(qin(:,3)) norm2delta(qin(:,4))];
    qbounds = [0 1; 0 15; 0 20; -2 2]';
elseif strcmp(modelID, 'ms_UG3_etaf_fixedNorm')
    if size(qin,2)~=4, disp('ERROR'); keyboard; end
    qout = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2f0(qin(:,3)) norm2delta(qin(:,4))];
    qbounds = [0 1; 0 15; 0 20; -2 2]';

   
end

