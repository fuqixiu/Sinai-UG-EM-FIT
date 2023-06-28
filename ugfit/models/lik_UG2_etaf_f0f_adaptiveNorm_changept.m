function [fval, norm, V, ChoiceProb] = lik_UG2_etaf_f0f_adaptiveNorm_changept(offer, resp, fixed, free, changept, doprior, varargin)
   %v2 prioritizes early information (trials < 20)
   %v3 change point to determines prioritization


    if nargin > 6
        prior      = varargin{1};
        q = varargin{2};
    end

    mn = 1;
    mx = 9;    

    eta = fixed(1);
    f0 = fixed(2);

    n = length(offer);

    [alpha, beta, epsilon, delta] = deal(free{:});


    % norm update (RW)
    norm = RW(f0, epsilon, offer);
    
   
    % value function / choice prob
    for i = 1:n       
        
        % consider 3 steps
        CV(i) = FS(alpha, offer(i), norm(i+1));   % net current value (accept - reject)        
        
        ao = max(offer(i)-delta, mn);
        if FS(alpha, ao, norm(i+1)) > 0          % if accept(now) & accept(next)
            aFV(i) = eta * FS(alpha, ao, norm(i+1)) + eta^2 * max(FS(alpha, max(ao-delta, mn), norm(i+1)), 0);
        else                                                                 % if accept & reject
            aFV(i) = eta^2 * max(FS(alpha, max(ao+delta, mn), norm(i+1)), 0);
        end
        
        ro = max(offer(i)+delta, mn);
        if FS(alpha, ro, norm(i+1)) > 0        % reject & accept
            rFV(i) = eta * FS(alpha, ro, norm(i+1)) + eta^2 * max(FS(alpha, max(ro-delta, mn), norm(i+1)), 0);
        else                                                                 % reject & reject
            rFV(i) = eta^2 * max(FS(alpha, max(ro+delta, mn), norm(i+1)), 0);
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


    log_probs = log(ChoiceProb);
    % V2: proportional discount of early/late trials
    early = 1:changept; %1:(n/2);
    late = (changept+1):n;
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