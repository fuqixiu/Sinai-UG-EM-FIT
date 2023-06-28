function [fval, norm, V, ChoiceProb] = lik_UG2_etaf_f0f_adaptiveNorm_blocked(offer_mat, resp_mat, fixed, free, doprior, varargin)
   if nargin > 5
        prior      = varargin{1};
        q = varargin{2};
    end

    mn = 1;
    mx = 9;    

    eta = fixed(1);
    f0 = fixed(2);

    n = size(offer_mat,1);
    b = size(offer_mat,2);

    resp = [resp_mat(:,1);  resp_mat(:,2)];

    [alpha, beta, epsilon, delta] = deal(free{:});

    
    V = zeros(n,b);
    % value function / choice prob
    for j = 1:b
        offer = offer_mat(:,j);
        % norm update (RW)
        norm = RW(f0, epsilon, offer);

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
                    
            V(i,j) = CV(i) + (aFV(i) - rFV(i));       
        end    
    end
        % calculate probability of accepting offer:                                                                                             
    prob_mat     = 1 ./ ( 1 + exp(-beta.*V));  
    
    prob = [prob_mat(:,1); prob_mat(:,2)];

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