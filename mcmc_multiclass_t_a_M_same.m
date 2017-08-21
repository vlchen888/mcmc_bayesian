function cont = mcmc_multiclass_t_a_M_same(params)
    % MCMC method, learning t, a, M
    % There is one copy of the hyperparameters shared between k classes
    % The pCN updates for xi are loop-based
    
    data = params('data');
    num_iterations = params('num_iterations');
    label_data = params('label_data');
    gamma = params('gamma');
    
    k = params('k');
    
    B_init = params('B_init');
    B_all = zeros(k, num_iterations);
    B_all(:,1) = B_init;
    
    init_tau = params('init_tau');
    init_alpha = params('init_alpha');
    init_M = params('init_M');
    
    tau_epsilon = params('tau_epsilon');
    alpha_epsilon = params('alpha_epsilon');
    M_max_jump = params('M_max_jump');
    
    tau_min = params('tau_min');
    tau_max = params('tau_max');
    alpha_min = params('alpha_min');
    alpha_max = params('alpha_max');
    M_min = params('M_min');
    M_max = params('M_max'); 
    
    tic;
    if params('laplacian') == "self tuning"
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == "un"
        p = params('p');
        q = params('q');
        l = params('l');
        L = compute_laplacian_standard(data, p, q, l);
    end
    
    fprintf('Laplacian computation complete: %.4f\n', toc);
    
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    if params('remove-zero-eig')
        lambda(1) = Inf;
    end
    
    [num_data, ~] = size(data);
    phi = phi(:, 1:M_max);
    lambda = lambda(1:M_max);
    tau_all = zeros(1, num_iterations);
    alpha_all = zeros(1, num_iterations);
    M_all = zeros(1, num_iterations);
    tau_all(1) = init_tau;
    alpha_all(1) = init_alpha;
    M_all(1) = init_M;
    
    xi_curr = zeros(M_max, k);
    u_avg = zeros(num_data, k);
    sign_avg = zeros(num_data, k);
    
        
    %%%%% Acceptance probabilities %%%%%
    xi_accept = zeros(k, num_iterations);
    tau_accept = zeros(1, num_iterations);
    alpha_accept = zeros(1, num_iterations);
    M_accept = zeros(1, num_iterations);
    
    uj_avg = zeros(M_max, k);
    
    %%%%% Store correct percent %%%%%
    correct_p = zeros(1, num_iterations);
    tic;    
    for i=1:num_iterations-1
        tau_curr = tau_all(i);
        alpha_curr = alpha_all(i);
        M_curr = M_all(i);
        u_curr = compute_T(xi_curr, tau_curr, alpha_curr, M_curr, lambda, phi);
        sign_curr = compute_S_multiclass(u_curr, k);
        if i >= params('burn_in')
            u_avg = (u_avg * (i-params('burn_in')) + u_curr)/(i-params('burn_in') + 1);
            sign_avg = (sign_avg * (i-params('burn_in')) + sign_curr)/(i-params('burn_in') + 1);
            
            uj_curr = (lambda + tau_curr^2).^(-alpha_curr/2) .* xi_curr / norm_const(lambda,tau_curr,alpha_curr);
            uj_curr(M_curr+1:end,:) = 0;
            uj_avg = (uj_avg * (i-params('burn_in')) + uj_curr)/(i-params('burn_in') + 1);
            if mod(i, params('movie_often')) == 0
                %p = count_correct_multiclass(compute_S_multiclass(u_avg, k), params('label_data'), params('truth'));            
                q = count_correct_multiclass(compute_S_multiclass(sign_avg, k), params('label_data'), params('truth'));

                correct_p(i) = q;
            end
        end
        
        %% Make figures
        if params('movie') && i >= params('burn_in') && mod(i, params('movie_often')) == 0
            fprintf('Sample number: %d, Time elapsed: %.2f\n', i, toc);
            %fprintf('Classification accuracy with S(E(u)): %.4f\n', p);
            fprintf('Classification accuracy with S(E(S(u))): %.4f\n', correct_p(i));
            fprintf('\txi accept acceptance probability: %.4f\n', mean(xi_accept(:,1:i),2));
            fprintf('\ttau accept acceptance probability: %.4f\n', mean(tau_accept(1:i)));
            fprintf('\talpha accept acceptance probability: %.4f\n', mean(alpha_accept(1:i)));
            fprintf('\tM accept acceptance probability: %.4f\n', mean(M_accept(1:i)));
            fprintf('tau running average: %.4f\n', mean(tau_all(1:i)));
            fprintf('alpha running average: %.4f\n', mean(alpha_all(1:i)));
            
            figure(1)
            set(gcf, 'Position', [0, 500, 500, 400])
            mnist_heatmap(compute_S_multiclass(sign_avg, k), params('truth'), params('digs'), "Confusion matrix, S(E(S(u)))");
            
            figure(2)
            set(gcf, 'Position', [500, 500, 700, 200])
            subplot(131)
            plot(u_avg)
            title('Average u')
            
            subplot(132)
            plot(uj_avg)
            title('Average u_j')
            
            subplot(133)
            plot(correct_p(correct_p~=0))
            title('Classification accuracy');
            xlabel(sprintf('Sample number (per %d)', params('movie_often')));
            
            figure(3)
            set(gcf, 'Position', [0, 0, 800, 300])
            subplot(121)
            plot(movmean(xi_accept(:,1:i)',[i 0]))
            title('\xi acceptance probability')
            
            subplot(122)
            plot(B_all(:,1:i)')
            title('\beta traces')
            
            figure(4)
            hyp_count = 0;
            fig_count = 1;
            if tau_epsilon > 0
                hyp_count = hyp_count + 1;
            end
            
            if alpha_epsilon > 0
                hyp_count = hyp_count + 1;
            end
            
            if params('M_max_jump') > 0
                hyp_count = hyp_count + 1;
            end
            
            set(gcf, 'Position', [500, 0, 700, 200*hyp_count])
            
            if tau_epsilon > 0
                subplot(hyp_count,2,fig_count)
                plot(1:i,tau_all(1:i),1:i,movmean(tau_all(1:i), [i 0]));
                legend('\tau trace', '\tau running average');

                subplot(hyp_count,2,fig_count+1)
                plot(movmean(tau_accept(1:i),[i 0]))
                title('\tau acceptance probability')
                
                fig_count = fig_count + 2;
            end
            
            if alpha_epsilon > 0
                subplot(hyp_count,2,fig_count)
                plot(1:i,alpha_all(1:i),1:i,movmean(alpha_all(1:i), [i 0]));
                legend('\alpha trace', '\alpha running average');

                subplot(hyp_count,2,fig_count+1)
                plot(movmean(alpha_accept(1:i),[i 0]))
                title('\alpha acceptance probability')
                
                fig_count = fig_count + 2;
            end
            
            if params('M_max_jump') > 0
                subplot(hyp_count,2,fig_count)
                plot(1:i, M_all(1:i), 1:i, movmean(M_all(1:i), [i 0]));
                legend('M trace', 'M running average');

                subplot(hyp_count,2,fig_count+1)
                plot(movmean(M_accept(1:i), [i 0]));
                title('M acceptance probability')
                
                %fig_count = fig_count + 2;
            end
            
            drawnow
        end
        %% xi proposal
        for j=1:k
            xi_hat = sqrt(1-B_all(j,i)^2)*xi_curr(:,j) + B_all(j,i)*normrnd(0,1,M_max,1);
            xi_new = xi_curr;
            xi_new(:, j) = xi_hat;
            u_old = compute_T(xi_curr, tau_curr, alpha_curr, M_curr, lambda, phi);
            u_new = compute_T(xi_new, tau_curr, alpha_curr, M_curr, lambda, phi);
            log_xi = compute_phi(gamma, label_data, u_old, k) - ...
                compute_phi(gamma, label_data, u_new, k);
            if rand(1) < exp(log_xi)
                xi_curr(:,j) = xi_hat;
                xi_accept(j,i+1) = 1;
            else
                xi_accept(j,i+1) = 0;
            end
        end
        
        %% tau proposal
        if tau_epsilon > 0
            tau_hat = tau_all(i) + tau_epsilon * normrnd(0, 1);
            if(tau_hat > tau_max) || (tau_hat < tau_min)
                tau_all(i+1) = tau_all(i);
                tau_accept(i+1) = 0;
            else
                u_old = compute_T(xi_curr, tau_curr, alpha_curr, M_curr, lambda, phi);
                u_hat = compute_T(xi_curr, tau_hat, alpha_curr, M_curr, lambda, phi);
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

        %% alpha proposal
        if alpha_epsilon > 0
            alpha_hat = alpha_all(i) + alpha_epsilon * normrnd(0, 1);
            if(alpha_hat > alpha_max) || (alpha_hat < alpha_min)
                alpha_all(i+1) = alpha_all(i);
                alpha_accept(i+1) = 0;
            else
                u_old = compute_T(xi_curr, tau_curr, alpha_curr, M_curr, lambda, phi);
                u_hat = compute_T(xi_curr, tau_curr, alpha_hat, M_curr, lambda, phi);
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
        
        %% M proposal
        if M_max_jump > 0
            M_hat = M_all(i) + compute_rndjump_M(M_max_jump);
            if (M_hat > M_max) || (M_hat < M_min)
                M_all(i+1) = M_all(i);
                M_accept(i+1) = 0;
            else
                u_old = compute_T(xi_curr, tau_curr, alpha_curr, M_curr, lambda, phi);
                u_hat = compute_T(xi_curr, tau_curr, alpha_curr, M_hat, lambda, phi);
                log_M = compute_phi(gamma, label_data, u_old, k) - compute_phi(gamma, label_data, u_hat, k);
                if rand(1) < exp(log_M)
                    M_all(i+1) = M_hat;
                    M_accept(i+1) = 1;
                else
                    M_all(i+1) = M_all(i);
                    M_accept(i+1) = 0;
                end
            end
        else
            M_all(i+1) = M_all(i);
            M_accept(i+1) = 0;
        end
        
        %% Updates to B
        if params('adaptive') && i <= params('B_burn_in') && mod(i, params('B_update_period')) == 0
            for j=1:k
                accept_mov_avg = mean(xi_accept(j, i-params('B_update_period')+1:i));
                %accept_tot_avg = mean(xi_accept(j, params('burn_in'):i));
                B_all(j,i+1) = min(B_all(j, i)*(1 + accept_mov_avg - params('B_target_p')), 1);
            end
        else
            B_all(:,i+1) = B_all(:,i);
        end
    end
    
    correct_p = correct_p(correct_p~=0);
    cont = containers.Map;
    cont('correct_p') = correct_p;
    %cont('M_all') = M_all;
    %cont('tau_all') = tau_all;
    %cont('alpha_all') = alpha_all;
end
%% Likelihood
function l = compute_phi(gamma, label_data, u, k)
    diff = abs(compute_S_multiclass(u, k) - label_data)/sqrt(2);
    diff = diff(sum(label_data, 2) ~= 0, :);
    l = sum(sum(diff))/(2*gamma^2);
end
%% Map from parameters to u
function T = compute_T(xi, tau, alpha, M, lambda, phi)
    T = phi(:,1:M)* ( (lambda(1:M) + tau^2).^(-alpha/2) .* xi(1:M,:)) / norm_const(lambda, tau, alpha);
end

%% M random jump
function j = compute_rndjump_M(k)
    prob_arr = zeros(2*k+1,1);
    for i = -k:k
        prob_arr(k+1+i) = 1/(1+abs(i));
    end
    prob_arr = prob_arr / sum(prob_arr);
    r = rand(1);
    for i=-k:k
        r = r - prob_arr(k+1+i);
        if r <= 0
            j = i;
            return
        end
    end        
end

%% Normalization const
function const = norm_const(lambda, tau, alpha)
    const = sqrt(sum((lambda+tau^2).^-alpha))/sqrt(length(lambda));
end