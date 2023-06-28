function [fval, norm, V, ChoiceProb] = lik_UG3_etaf_fixedNorm(offer, resp, fixed, free, doprior, varargin)
    if nargin > 5
        prior      = varargin{1};
        q = varargin{2};
    end
    mn = 1;    

    n = length(offer);

    eta = fixed(1);

    [alpha, beta, f0, delta] = deal(free{:});

    % norm update (RW)
    norm = f0; %RW(f0, epsilon, offer);
    
   
    % value function / choice prob
    for i = 1:n               
        
        CV(i) = FS( alpha, offer(i), norm);   % net current value (accept - reject)        
        
        
        % FV_accept        
        
        ao = max(offer(i)-delta, mn);       % ao = expected offer1 when accept offer0
        
        if FS( alpha, ao, norm) > 0                          % a a            
                aao = max(ao-delta, mn);                
                if FS( alpha, aao, norm) > 0                 % a a a
                    aFV(i) = eta * FS( alpha, ao, norm) + eta^2 * FS( alpha, aao, norm) ...
                        + eta^3 * max(FS( alpha, max(aao-delta, mn), norm), 0);
                else                                                            % a a r
                    aFV(i) = eta * FS( alpha, ao, norm) + eta^2 * 0 ...
                        + eta^3 * max(FS( alpha, max(aao+delta, mn), norm), 0);
                end
            
        else                                                                 % a r
                aro = max(ao+delta, mn);                
                if FS( alpha, aro, norm) > 0                 % a r a
                    aFV(i) = eta * 0 + eta^2 * FS( alpha, aro, norm) ...
                        + eta^3 * max(FS( alpha, max(aro-delta, mn), norm), 0);
                else                                                            % a r r
                    aFV(i) = eta * 0 + eta^2 * 0 ...
                        + eta^3 * max(FS( alpha, max(aro+delta, mn), norm), 0);
                end
        end
        
        
        % FV_reject
        
        ro = max(offer(i)+delta, mn);       % ro = expected offer1 when reject offer0
        
        if FS( alpha, ro, norm) > 0                          % r a            
                rao = max(ro-delta, mn);                
                if FS( alpha, rao, norm) > 0                 % r a a
                    rFV(i) = eta * FS( alpha, ro, norm) + eta^2 * FS( alpha, rao, norm) ...
                        + eta^3 * max(FS( alpha, max(rao-delta, mn), norm), 0);
                else                                                            % r a r
                    rFV(i) = eta * FS( alpha, ro, norm) + eta^2 * 0 ...
                        + eta^3 * max(FS( alpha, max(rao+delta, mn), norm), 0);
                end
            
        else                                                                 % r r
                rro = max(ao+delta, mn);                
                if FS( alpha, rro, norm) > 0                 % r r a
                    rFV(i) = eta * 0 + eta^2 * FS( alpha, rro, norm) ...
                        + eta^3 * max(FS( alpha, max(rro-delta, mn), norm), 0);
                else                                                            % r r r
                    rFV(i) = eta * 0 + eta^2 * 0 ...
                        + eta^3 * max(FS( alpha, max(rro+delta, mn), norm), 0);
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

    nll =-nansum(log(ChoiceProb));  % the thing to minimize                      
    
    if doprior == 0                % NLL fit
       fval = nll;
    elseif doprior == 1             % EM-fit:   P(Choices | h) * P(h | O) should be maximised, therefore same as minimizing it with negative sign   
       fval = -(-nll + prior.logpdf(q));
    end
    
    if sum(isnan(ChoiceProb))>0, disp('ERROR'); keyboard; return; end   % error check 

end