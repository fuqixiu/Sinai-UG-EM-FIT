function [fval, norm, V, ChoiceProb] = lik_UG2_etaf_fixedNorm(offer, resp, fixed, free, doprior, varargin)
   if nargin > 5
        prior      = varargin{1};
        q = varargin{2};
    end

    mn = 1;
    mx = 9;    

    n = length(offer);

    eta = fixed(1);

    n = length(offer);

    [alpha, beta, f0, delta] = deal(free{:});


    % norm update (RW)
    norm = f0;
   
    % value function / choice prob
    for i = 1:n       
        
        % consider 3 steps
        CV(i) = FS(alpha, offer(i), norm);   % net current value (accept - reject)        
        
        ao = max(offer(i)-delta, mn);
        if FS(alpha, ao, norm) > 0          % if accept(now) & accept(next)
            aFV(i) = eta * FS(alpha, ao, norm) + eta^2 * max(FS(alpha, max(ao-delta, mn), norm), 0);
        else                                                                 % if accept & reject
            aFV(i) = eta^2 * max(FS(alpha, max(ao+delta, mn), norm), 0);
        end
        
        ro = max(offer(i)+delta, mn);
        if FS(alpha, ro, norm) > 0        % reject & accept
            rFV(i) = eta * FS(alpha, ro, norm) + eta^2 * max(FS(alpha, max(ro-delta, mn), norm), 0);
        else                                                                 % reject & reject
            rFV(i) = eta^2 * max(FS(alpha, max(ro+delta, mn), norm), 0);
        end
                
        V(i) = CV(i) + (aFV(i) - rFV(i));       
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