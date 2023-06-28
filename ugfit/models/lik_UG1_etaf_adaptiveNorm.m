function [fval, norm, V, ChoiceProb] = lik_UG1_etaf_adaptiveNorm(offer, resp, fixed, free, doprior, varargin)
    mn = 1;  

    if nargin > 5
        prior      = varargin{1};
        q = varargin{2};
    end
    
    n = length(offer);

    eta = fixed(1);
    
   [alpha, beta, f0, epsilon, delta] = deal(free{:});
   

    % norm update (RW)
    norm = RW(f0, epsilon, offer);
    
   
    % value function / choice prob
    for i = 1:n       
        
        CV(i) =FS(alpha, offer(i), norm(i+1));
        FVa(i) = max(0, FS(alpha,  max(offer(i)-delta, mn), norm(i+1)));
        FVr(i) = max(0, FS(alpha,  max(offer(i)+delta, mn), norm(i+1)));
        V(i) = CV(i) + eta * ( FVa(i) - FVr(i) );
        
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