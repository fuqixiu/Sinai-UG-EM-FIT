% Adapted by BRK Shevlin April 2023
% Adapted by AD August 2023

function [fval, norm, V, ChoiceProb] = lik_UG0_adaptiveNorm_fixedalpha(offer, resp, ~, free, doprior, varargin)

    if nargin > 5
        prior      = varargin{1};
        q = varargin{2};
    end

    n = length(offer);

    %[alpha, beta, f0, epsilon, ~] = deal(free{:});
    [beta, f0, epsilon] = deal(free{:});
    alpha = 0.8;

    % norm update (RW)
    norm = RW(f0, epsilon, offer);  

    % value function / choice prob
    for i = 1:n        
        V(i) = FS(alpha, offer(i), norm(i+1));   % net current value (accept - reject)              
        
        %Like = Like - log(exp(beta*V(i)*(resp(i)))/(1+exp(beta*V(i)))); 
    end

    % calculate probability of accepting offer:                                                                                             
    prob     = 1 ./ ( 1 + exp(-beta.*V));  

    % find when offer was actually choosen:
    accept = find(resp == 1);
    reject = find(resp == 0);

    ChoiceProb(accept) = prob(accept);
    ChoiceProb(reject) = 1 - prob(reject);

    nll =-nansum(log(ChoiceProb));  % the thing to minimize                      
    
    if doprior == 0                % NLL fit
       fval = nll;
    elseif doprior == 1             % EM-fit:   P(Choices | h) * P(h | O) should be maximised, therefore same as minimizing it with negative sign   
       fval = -(-nll + prior.logpdf(q));
    end
    
    if sum(isnan(ChoiceProb))>0, disp('ERROR'); keyboard; return; end   % error check  

end