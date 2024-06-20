function [fval,fit] = mod_ms_UG0_adaptiveNorm_fixedalpha(behavData,q, doprior,dofit,varargin)
% runs standard 0 Step UG Model
% Created MK Wittmann, Oct 2018
% Adapted BRK Shevlin, April 2023
% Adapted AD August 2023
%
% INPUT:    - behavData: behavioural input file
% OUTPUT:   - fval and fitted variables
% 
% 
%
%%
% -------------------------------------------------------------------------------------
% 1 ) Define free parameters
% -------------------------------------------------------------------------------------

if nargin > 4
    prior      = varargin{1};
end

[qt, bounds] = norm2par('ms_UG0_adaptiveNorm_fixedalpha',q); % transform parameters from gaussian space to model space

% Define free parameters and set unused ones to zero
beta = qt(1); % envy
if (beta<min(bounds(:,1)) || beta>max(bounds(:,1))), fval=10000000; return; end

f0 = qt(2); % initial norm
if (f0<min(bounds(:,2)) || f0>max(bounds(:,2))), fval=10000000; return; end

epsilon = qt(3); % norm adaptation rate
if (epsilon<min(bounds(:,3)) || epsilon>max(bounds(:,3))), fval=10000000; return; end

fixed = 0.8;
free = {beta f0 epsilon};
% -------------------------------------------------------------------------------------
% 2-4) Middle code is specific the the model
% -------------------------------------------------------------------------------------
if doprior == 1
    [fval,norm,V,ChoiceProb] = lik_UG0_adaptiveNorm_fixedalpha(behavData.offer, behavData.choice,fixed,free,doprior,prior,q);
else
    [fval,norm,V,ChoiceProb] = lik_UG0_adaptiveNorm_fixedalpha(behavData.offer, behavData.choice,fixed,free,doprior);

end
% -------------------------------------------------------------------------------------
% 5) Calculate additional Parameters and save: 
% -------------------------------------------------------------------------------------

if dofit ==1

   fit         = struct;
   fit.xnames  = {'beta';'f0';'epsilon'};
   
   fit.mat    = [ChoiceProb norm V];
   fit.names  = {'ChoiceProb' ; 'norm'; 'V'};
end




end

