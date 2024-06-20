% Adapted by BRK Shevlin April 2023

function [fval, norm, V, ChoiceProb] = lik_UG0_Bayes_10(offer, resp, fixed, free, doprior, varargin)

    if nargin > 5
        prior      = varargin{1};
        q = varargin{2};
    end

    n = length(offer);

    [alpha, beta] = deal(free{:});

    f0 = fixed(1);
    k = fixed(2);
    norm = zeros(1, n+1);
    norm(1) = f0;

    % norm update (Bayes)
    for i = 1:n
        k = k+1;
        norm(1+i) = (k-1) / k * norm(i) + 1 / k * offer(i);
    end 

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
       if isinf(fval)
           fval = 10000000;
           return
       end
    end
    
    if sum(isnan(ChoiceProb))>0, disp('ERROR'); keyboard; return; end   % error check  

end