function [u_all] = mcmc_multiclass_t_a_same(params)
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
    u_all = zeros(num_data, k, num_iterations);
        
    %%%%% Acceptance probabilities %%%%%
    xi_accept = zeros(k, num_iterations);
    tau_accept = zeros(1, num_iterations);
    alpha_accept = zeros(1, num_iterations);
    
    tic;    
    for i=1:num_iterations-1
        tau_curr = tau_all(i);
        alpha_curr = alpha_all(i);
        %signs_all(:,:,i) = compute_S(compute_T(xi_curr, tau, alpha, lambda, phi), k);
        u_all(:,:,i) = compute_T(xi_curr, tau_curr, alpha_curr, lambda, phi);
        
        if i >= params('burn_in') && mod(i, 2500) == 0
            
            avg_label = zeros(num_data, k);
            for ii = params('burn_in') : i
                avg_label = avg_label + compute_S_multiclass(squeeze(u_all(:,:,ii)), k);
            end
            %avg_label = mean(u_all(:,:,params('burn_in'):i),3);
            curr_label = compute_S_multiclass(avg_label, k);
            p = count_correct_multiclass(curr_label, params('label_data'), params('truth'));
            fprintf('Sample number: %d, Time elapsed: %.2f\n', i, toc);
            fprintf('Classification accuracy: %.4f\n', p);
            fprintf('\txi accept acceptance probability: %.4f\n', mean(xi_accept(:,1:i),2));
            fprintf('\ttau accept acceptance probability: %.4f\n', mean(tau_accept(1:i),2));
            fprintf('\talpha accept acceptance probability: %.4f\n', mean(alpha_accept(1:i),2));
            fprintf('tau running average: %.4f\n', mean(tau_all(1:i)));
            fprintf('alpha running average: %.4f\n', mean(alpha_all(1:i)));
            figure(1)
            mnist_heatmap(curr_label, params('truth'), params('digs'));
            
            figure(2)
            subplot(211)
            plot(u_all(:,:,i))
            subplot(212)
            plot((lambda + tau_curr^2).^(-alpha_curr/2) .* xi_curr)
            
            figure(3)
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
    diff = diff( sum(label_data, 2) ~= 0 , : );
    l = sum(sum(diff))/(2*gamma^2);
end

function T = compute_T(xi, tau, alpha, lambda, phi)
    T = phi* ( (lambda + tau^2).^(-alpha/2) .* xi);
end