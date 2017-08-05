function [u_avg] = mcmc_multiclass_t_a(params)
    data = params('data');
    num_iterations = params('num_iterations');
    label_data = params('label_data');
    p = params('p');
    q = params('q');
    l = params('l');
    gamma = params('gamma');
    B = params('B');
    k = params('k');
    
    
    init_tau = params('init_tau');
    init_alpha = params('init_alpha');
    
    tau_epsilon = params('tau_epsilon');
    alpha_epsilon = params('alpha_epsilon');
    
    tau_min = params('tau_min');
    tau_max = params('tau_max');
    alpha_min = params('alpha_min');
    alpha_max = params('alpha_max');
    
    tic;
    if params('laplacian') == "self tuning"
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == "un"
        L = compute_laplacian_standard(data, p, q, l);
    end
    
    fprintf('Laplacian computation complete: %.4f\n', toc);
    
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    
    [num_data, ~] = size(data);
    M = 50;
    phi = phi(:, 1:M);
    lambda = lambda(1:M);
    
    xi_curr = zeros(M, k);
    u_avg = zeros(num_data, k);
    sign_avg = zeros(num_data, k);
    taus_curr = zeros(1, k);
    alphas_curr = zeros(1, k);
    
    taus_curr(1,:) = init_tau;
    alphas_curr(1,:) = init_alpha;
    
    taus_all = zeros(k, num_iterations);
    alphas_all = zeros(k, num_iterations);
    
    %%%%% Acceptance probabilities %%%%%
    xi_accept = zeros(k, num_iterations);
    tau_accept = zeros(k, num_iterations);
    alpha_accept = zeros(k, num_iterations);
    
    tic;    
    for i=1:num_iterations-1
        taus_all(:, i) = taus_curr;
        alphas_all(:, i) = alphas_curr;
        u_curr = compute_T(xi_curr, taus_curr, alphas_curr, lambda, phi);
        sign_curr = compute_S_multiclass(u_curr, k);
        
        if i >= params('burn_in')
            u_avg = (u_avg * (i-params('burn_in')) + u_curr)/(i-params('burn_in') + 1);
            sign_avg = (sign_avg * (i-params('burn_in')) + sign_curr)/(i-params('burn_in') + 1);
        end
        
        if i >= params('burn_in') && mod(i, 2500) == 0
            
            p = count_correct_multiclass(compute_S_multiclass(u_avg, k), params('label_data'), params('truth'));
            q = count_correct_multiclass(compute_S_multiclass(sign_avg, k), params('label_data'), params('truth'));
            fprintf('Sample number: %d, Time elapsed: %.2f\n', i, toc);
            fprintf('Classification accuracy with S(E(u)): %.4f\n', p);
            fprintf('Classification accuracy with S(E(S(u))): %.4f\n', q);
            fprintf('\txi accept acceptance probability: %.4f\n', mean(xi_accept(:,1:i),2));
            fprintf('\ttau accept acceptance probability: %.4f\n', mean(tau_accept(:,1:i),2));
            fprintf('\talpha accept acceptance probability: %.4f\n', mean(alpha_accept(:,1:i),2));
            figure(1)
            mnist_heatmap(compute_S_multiclass(sign_avg, k), params('truth'), params('digs'), "S(E(S(u))) classification");
            
            figure(3)
            set(gcf, 'Position', [0, 0, 500, 300])
            subplot(311)
            plot(u_curr)
            title('Current u')
            subplot(312)
            plot((lambda + taus_curr.^2).^(-alphas_curr/2) .* xi_curr)
            title('Current u_j')
            subplot(313)
            plot(u_avg);
            title('Average u')
            
            figure(4)
            set(gcf, 'Position', [500, 0, 500, 300])
            for kk = 1:k
                subplot(k,1,kk)
                plot(movmean(xi_accept(kk,1:i),[i 0]))
                if kk == 1
                    title('\xi acceptance probability')
                end
            end
            
            figure(5)
            set(gcf, 'Position', [0, 0, 500, 300])
            subplot(211)
            plot(taus_all(:,1:i))
            plot(1:i,taus_all(:,1:i)',1:i,movmean(taus_all(:,1:i)',[i 0]))
            xlabel('\tau trace')
            subplot(212)
            plot(1:i,alphas_all(:,1:i)',1:i,movmean(alphas_all(:,1:i)',[i 0]))
            xlabel('\alpha trace')
            
            drawnow
        end
        
        for j=1:k
            %xi proposal
            xi_hat = sqrt(1-B^2)*xi_curr(:,j) + B*normrnd(0,1,M,1);
            xi_new = xi_curr;
            xi_new(:, j) = xi_hat;
            u_old = compute_T(xi_curr, taus_curr, alphas_curr, lambda, phi);
            u_new = compute_T(xi_new, taus_curr, alphas_curr, lambda, phi);
            log_xi = compute_phi(gamma, label_data, u_old, k) - ...
                compute_phi(gamma, label_data, u_new, k);
            if rand(1) < exp(log_xi)
                xi_curr(:,j) = xi_hat;
                xi_accept(j,i+1) = 1;
            else
                xi_accept(j,i+1) = 0;
            end
            
            %tau proposal
            tau_hat = taus_curr(j) + tau_epsilon * normrnd(0, 1);
            if(tau_hat > tau_max) || (tau_hat < tau_min)
                tau_accept(j,i+1) = 0;
            else
                new_taus = taus_curr;
                new_taus(1, j) = tau_hat;
                u_old = compute_T(xi_curr, taus_curr, alphas_curr, lambda, phi);
                u_hat = compute_T(xi_curr, new_taus, alphas_curr, lambda, phi);
                log_tau = compute_phi(gamma, label_data, u_old, k) - compute_phi(gamma, label_data, u_hat, k);
                
                if rand(1) < exp(log_tau)
                    taus_curr(1, j) = tau_hat;
                    tau_accept(j,i+1) = 1;
                else
                    tau_accept(j,i+1) = 0;
                end
            end
            
            %alpha proposal
            alpha_hat = alphas_curr(j) + alpha_epsilon * normrnd(0, 1);
            if(alpha_hat > alpha_max) || (alpha_hat < alpha_min)
                alpha_accept(j,i+1) = 0;
            else
                new_alphas = alphas_curr;
                new_alphas(1, j) = alpha_hat;
                u_old = compute_T(xi_curr, taus_curr, alphas_curr, lambda, phi);
                u_hat = compute_T(xi_curr, taus_curr, new_alphas, lambda, phi);
                log_alpha = compute_phi(gamma, label_data, u_old, k) - compute_phi(gamma, label_data, u_hat, k);
                
                if rand(1) < exp(log_alpha)
                    alphas_curr(1, j) = alpha_hat;
                    alpha_accept(j,i+1) = 1;
                else
                    alpha_accept(j,i+1) = 0;
                end
            end
        end
    end
    
end

function l = compute_phi(gamma, label_data, u, k)
    diff = abs(compute_S_multiclass(u, k) - label_data)/sqrt(2);
    diff = diff(sum(label_data, 2) ~= 0, : );
    l = sum(sum(diff))/(2*gamma^2);
end

function T = compute_T(xi, taus, alphas, lambda, phi)
    T = phi* ( (lambda + taus.^2).^(-alphas/2) .* xi);
end