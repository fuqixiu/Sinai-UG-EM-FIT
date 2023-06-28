function [fval, norm, V, ChoiceProb] = lik_UG3_etaf_f0f_adaptiveNorm_v2(offer, resp, fixed, free, doprior, varargin)
    if nargin > 5
        prior      = varargin{1};
        q = varargin{2};
    end
    mn = 1;    

    n = length(offer);

    eta = fixed(1);
    f0 = fixed(2);

    [alpha, beta, epsilon, delta] = deal(free{:});

    % norm update (RW)
    norm = RW(f0, epsilon, offer);
    
   
    % value function / choice prob
    for i = 1:n               
        
        CV(i) = FS( alpha, offer(i), norm(i+1));   % net current value (accept - reject)        
        
        
        % FV_accept        
        
        ao = max(offer(i)-delta, mn);       % ao = expected offer1 when accept offer0
        
        if FS( alpha, ao, norm(i+1)) > 0                          % a a            
                aao = max(ao-delta, mn);                
                if FS( alpha, aao, norm(i+1)) > 0                 % a a a
                    aFV(i) = eta * FS( alpha, ao, norm(i+1)) + eta^2 * FS( alpha, aao, norm(i+1)) ...
                        + eta^3 * max(FS( alpha, max(aao-delta, mn), norm(i+1)), 0);
                else                                                            % a a r
                    aFV(i) = eta * FS( alpha, ao, norm(i+1)) + eta^2 * 0 ...
                        + eta^3 * max(FS( alpha, max(aao+delta, mn), norm(i+1)), 0);
                end
            
        else                                                                 % a r
                aro = max(ao+delta, mn);                
                if FS( alpha, aro, norm(i+1)) > 0                 % a r a
                    aFV(i) = eta * 0 + eta^2 * FS( alpha, aro, norm(i+1)) ...
                        + eta^3 * max(FS( alpha, max(aro-delta, mn), norm(i+1)), 0);
                else                                                            % a r r
                    aFV(i) = eta * 0 + eta^2 * 0 ...
                        + eta^3 * max(FS( alpha, max(aro+delta, mn), norm(i+1)), 0);
                end
        end
        
        
        % FV_reject
        
        ro = max(offer(i)+delta, mn);       % ro = expected offer1 when reject offer0
        
        if FS( alpha, ro, norm(i+1)) > 0                          % r a            
                rao = max(ro-delta, mn);                
                if FS( alpha, rao, norm(i+1)) > 0                 % r a a
                    rFV(i) = eta * FS( alpha, ro, norm(i+1)) + eta^2 * FS( alpha, rao, norm(i+1)) ...
                        + eta^3 * max(FS( alpha, max(rao-delta, mn), norm(i+1)), 0);
                else                                                            % r a r
                    rFV(i) = eta * FS( alpha, ro, norm(i+1)) + eta^2 * 0 ...
                        + eta^3 * max(FS( alpha, max(rao+delta, mn), norm(i+1)), 0);
                end
            
        else                                                                 % r r
                rro = max(ao+delta, mn);                
                if FS( alpha, rro, norm(i+1)) > 0                 % r r a
                    rFV(i) = eta * 0 + eta^2 * FS( alpha, rro, norm(i+1)) ...
                        + eta^3 * max(FS( alpha, max(rro-delta, mn), norm(i+1)), 0);
                else                                                            % r r r
                    rFV(i) = eta * 0 + eta^2 * 0 ...
                        + eta^3 * max(FS( alpha, max(rro+delta, mn), norm(i+1)), 0);
                end
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