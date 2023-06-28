% Adapted by BRK Shevlin April 2023
% May 2023: v2 weights early trials > late trials

function [fval, norm, V, ChoiceProb] = lik_UG0_f0f_adaptiveNorm_v2(offer, resp, fixed, free, doprior, varargin)

    if nargin > 5
        prior      = varargin{1};
        q = varargin{2};
    end

    n = length(offer);

    [alpha, beta, epsilon] = deal(free{:});

    f0 = fixed(1);

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

    % V2: proportional discount of early/late trials
    log_probs = log(ChoiceProb);
    early = 1:(n/2);
    late = (n/2)+1:n;
    log_probs(early) = log_probs(early) * 4;
    log_probs(late) = log_probs(late);
    nll =-nansum(log_probs);  % the thing to minimize 

    
    if doprior == 0                % NLL fit
       fval = nll;
    elseif doprior == 1             % EM-fit:   P(Choices | h) * P(h | O) should be maximised, therefore same as minimizing it with negative sign   
       fval = -(-nll + prior.logpdf(q));
       if isinf(fval)
           fval = 10000000;
           return
       end
    end
    
    if sum(isnan(ChoiceProb))>0, disp('ERROR'); keyboard; return; end   % error check  

end