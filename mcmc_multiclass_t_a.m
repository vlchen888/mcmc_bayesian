function [u_all] = mcmc_multiclass_t_a(params)
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
    u_all = zeros(num_data, k, num_iterations);
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
        %signs_all(:,:,i) = compute_S(compute_T(xi_curr, tau, alpha, lambda, phi), k);
        taus_all(:, i) = taus_curr;
        alphas_all(:, i) = alphas_curr;
        u_all(:,:,i) = compute_T(xi_curr, taus_curr, alphas_curr, lambda, phi);
        
        if mod(i, 2500) == 0
            avg_label = mean(u_all, 3);
            curr_label = compute_S(avg_label, k);
            p = count_correct_multiclass(curr_label, params('label_data'), params('truth'));
            fprintf('Sample number: %d, Time elapsed: %.2f\n', i, toc);
            fprintf('Classification accuracy: %.4f\n', p);
            fprintf('\txi accept acceptance probability: %.4f\n', mean(xi_accept(:,1:i),2));
            fprintf('\ttau accept acceptance probability: %.4f\n', mean(tau_accept(:,1:i),2));
            fprintf('\talpha accept acceptance probability: %.4f\n', mean(alpha_accept(:,1:i),2));
            figure(1)
            mnist_heatmap(curr_label, params('truth'), params('digs'));
            
            figure(2)
            plot(u_all(:,:,i))
            
            figure(3)
            subplot(211)
            plot(taus_all(:,1:i)')
            xlabel('\tau s')
            subplot(212)
            plot(alphas_all(:,1:i)')
            xlabel('\alpha s')
            
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
    diff = abs(compute_S(u, k) - label_data)/sqrt(2);
    diff = diff(:, find (sum(label_data) ~= 0) );
    l = sum(sum(diff))/(2*gamma^2);
end

function S = compute_S(u, k)
    [~, I] = max(u,[],2);
    A = eye(k);
    S = A(:, I);
end

function T = compute_T(xi, taus, alphas, lambda, phi)
    T = phi* ( (lambda + taus.^2).^(-alphas/2) .* xi);
end