function [u_avg, sign_avg] = mcmc_multiclass_t_a_same(params)
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
    tau_all = zeros(1, num_iterations);
    alpha_all = zeros(1, num_iterations);
    tau_all(1) = init_tau;
    alpha_all(1) = init_alpha;
    
    xi_curr = zeros(M, k);
    u_avg = zeros(num_data, k);
    sign_avg = zeros(num_data, k);
    
        
    %%%%% Acceptance probabilities %%%%%
    xi_accept = zeros(k, num_iterations);
    tau_accept = zeros(1, num_iterations);
    alpha_accept = zeros(1, num_iterations);
    
    tic;    
    for i=1:num_iterations-1
        tau_curr = tau_all(i);
        alpha_curr = alpha_all(i);
        u_curr = phi * ((lambda + tau_curr^2).^(-alpha_curr/2) .* xi_curr);
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
            fprintf('\ttau accept acceptance probability: %.4f\n', mean(tau_accept(1:i),2));
            fprintf('\talpha accept acceptance probability: %.4f\n', mean(alpha_accept(1:i),2));
            fprintf('tau running average: %.4f\n', mean(tau_all(1:i)));
            fprintf('alpha running average: %.4f\n', mean(alpha_all(1:i)));
            
            figure(1)
            set(gcf, 'Position', [0, 500, 500, 400])
            mnist_heatmap(compute_S_multiclass(sign_avg, k), params('truth'), params('digs'), "S(E(S(u))) classification");
            
            figure(3)
            set(gcf, 'Position', [500, 500, 500, 700])
            subplot(411)
            plot(u_curr)
            title('Current u')
            subplot(412)
            plot((lambda + tau_curr^2).^(-alpha_curr/2) .* xi_curr)
            title('Current u_j')
            subplot(413)
            plot(u_avg);
            title('Average u')
            subplot(414)
            plot(sign_avg);
            title('Average S(u)')
            
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
            plot(tau_all(1:i))
            plot(1:i,tau_all(1:i),1:i,movmean(tau_all(1:i),[i 0]))
            xlabel('\tau trace')
            subplot(212)
            plot(1:i,alpha_all(1:i),1:i,movmean(alpha_all(1:i),[i 0]))
            xlabel('\alpha trace')
            
            drawnow
        end
        
        for j=1:k
            %xi proposal
            xi_hat = sqrt(1-B^2)*xi_curr(:,j) + B*normrnd(0,1,M,1);
            xi_new = xi_curr;
            xi_new(:, j) = xi_hat;
            u_old = compute_T(xi_curr, tau_curr, alpha_curr, lambda, phi);
            u_new = compute_T(xi_new, tau_curr, alpha_curr, lambda, phi);
            log_xi = compute_phi(gamma, label_data, u_old, k) - ...
                compute_phi(gamma, label_data, u_new, k);
            if rand(1) < exp(log_xi)
                xi_curr(:,j) = xi_hat;
                xi_accept(j,i+1) = 1;
            else
                xi_accept(j,i+1) = 0;
            end
        end
        
        %tau proposal
        if tau_epsilon > 0
            tau_hat = tau_all(i) + tau_epsilon * normrnd(0, 1);
            if(tau_hat > tau_max) || (tau_hat < tau_min)
                tau_all(i+1) = tau_all(i);
                tau_accept(i+1) = 0;
            else
                u_old = compute_T(xi_curr, tau_curr, alpha_curr, lambda, phi);
                u_hat = compute_T(xi_curr, tau_hat, alpha_curr, lambda, phi);
                log_tau = compute_phi(gamma, label_data, u_old, k) - compute_phi(gamma, label_data, u_hat, k);

                if rand(1) < exp(log_tau)
                    tau_all(i+1) = tau_hat;
                    tau_accept(i+1) = 1;
                else
                    tau_all(i+1) = tau_all(i);
                    tau_accept(i+1) = 0;
                end
            end
        else
            tau_all(i+1) = tau_all(i);
            tau_accept(i+1) = 0;
        end

        %alpha proposal
        if alpha_epsilon > 0
            alpha_hat = alpha_all(i) + alpha_epsilon * normrnd(0, 1);
            if(alpha_hat > alpha_max) || (alpha_hat < alpha_min)
                alpha_all(i+1) = alpha_all(i);
                alpha_accept(i+1) = 0;
            else
                u_old = compute_T(xi_curr, tau_curr, alpha_curr, lambda, phi);
                u_hat = compute_T(xi_curr, tau_curr, alpha_hat, lambda, phi);
                log_alpha = compute_phi(gamma, label_data, u_old, k) - compute_phi(gamma, label_data, u_hat, k);

                if rand(1) < exp(log_alpha)
                    alpha_all(i+1) = alpha_hat;
                    alpha_accept(i+1) = 1;
                else
                    alpha_all(i+1) = alpha_all(i);
                    alpha_accept(i+1) = 0;
                end
            end
        else
            alpha_all(i+1) = alpha_all(i);
            alpha_accept(i+1) = 0;
        end
    end
    
end

function l = compute_phi(gamma, label_data, u, k)
    diff = abs(compute_S_multiclass(u, k) - label_data)/sqrt(2);
    diff = diff(sum(label_data, 2) ~= 0, :);
    l = sum(sum(diff))/(2*gamma^2);
end

function T = compute_T(xi, tau, alpha, lambda, phi)
    T = phi* ( (lambda + tau^2).^(-alpha/2) .* xi);
end