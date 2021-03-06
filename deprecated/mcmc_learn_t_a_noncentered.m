function [tau_all, alpha_all, std, xi_accept, tau_accept, alpha_accept] =...
        mcmc_learn_t_a_noncentered(params)
    data = params('data');
    num_iterations = params('num_iterations');
    label_data = params('label_data');
    p = params('p');
    q = params('q');
    l = params('l');
    gamma = params('gamma');
    B = params('B');
    init_tau = params('init_tau');
    init_alpha = params('init_alpha');
    min_tau = params('min_tau');
    max_tau = params('max_tau');
    min_alpha = params('min_alpha');
    max_alpha = params('max_alpha');
    alpha_epsilon = params('alpha_epsilon');
    tau_epsilon = params('tau_epsilon');
    
    if params('laplacian') == "self tuning"
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == "un"
        L = compute_laplacian_standard(data, p, q, l);
    end
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    
    [num_data, ~] = size(data);
    M = 50;
    phi = phi(:, 1:M);
    lambda = lambda(1:M);
    
    xi_all = zeros(M, num_iterations);
    
    %%%%% Initialization from Fiedler Vector?? %%%%%
    %xi_all(2, 1) = (lambda(2)+init_tau^2)^(-init_alpha/2);
        
    tau_all = zeros(1, num_iterations);
    tau_all(1) = init_tau;
    
    alpha_all = zeros(1, num_iterations);
    alpha_all(1) = init_alpha;
    
    %%%%% Acceptance probabilities %%%%%
    xi_accept       = zeros(1, num_iterations);
    tau_accept      = zeros(1, num_iterations);
    alpha_accept    = zeros(1, num_iterations);
    
    %%%%% Store standard basis vectors %%%%%
    std = zeros(num_data, num_iterations);
    std(:, 1) = compute_T(xi_all(:, 1), init_tau, init_alpha, lambda, phi);
    
    uj_all = zeros(M, num_iterations);
        
    for i=1:num_iterations-1
        
        if mod(i,2500) == 0
            figure(1)
            u_avg = mean(std(:,1:i), 2);
            plot(u_avg)
            drawnow
            p = count_correct(u_avg, params('label_data'), params('truth'));
            fprintf('Sample number: %d\n', i);
            fprintf('Classification accuracy: %.4f\n', p);
            fprintf('\txi accept acceptance probability: %.4f\n', mean(xi_accept(1:i)));
        end
        
        %%%%% Propose new state for xi %%%%%
        curr_xi = xi_all(:,i);
        curr_tau = tau_all(i);
        curr_alpha = alpha_all(i);
        
        uj_all(:, i) = (lambda + curr_tau^2).^(-curr_alpha/2) .* curr_xi;
        
        std(:, i) = compute_T(curr_xi, curr_tau, curr_alpha, lambda, phi);
        
        x = compute_rand_xi(M);
        new_xi = (1-B^2)^0.5*curr_xi+B*x;
        u_curr = compute_T(curr_xi, curr_tau, curr_alpha, lambda, phi);
        u_new = compute_T(new_xi, curr_tau, curr_alpha, lambda, phi);
        log_xi_trans = compute_phi(gamma, label_data, u_curr)...
            - compute_phi(gamma, label_data, u_new);
        transition_xi = exp(log_xi_trans);
        
        %%%% Do transition %%%%
        if rand(1) < transition_xi
            xi_all(:, i+1) = new_xi;
            xi_accept(i+1) = 1;
        else
            xi_all(:, i+1) = xi_all(:, i);
            xi_accept(i+1) = 0;
        end
        
        curr_xi = xi_all(:,i+1);
                
        %%%%% Propose a new tau %%%%%
        new_tau = tau_all(i) + tau_epsilon * normrnd(0, 1);
        if new_tau < min_tau || new_tau > max_tau
            tau_all(i+1) = tau_all(i);
            tau_accept(i+1) = 0;
        else
            log_tau = compute_log_g(lambda, phi, curr_xi, new_tau, curr_alpha, gamma, label_data) ...
            - compute_log_g(lambda, phi, curr_xi, curr_tau, curr_alpha, gamma, label_data);
            transition_tau = exp(log_tau);
            if rand(1) < transition_tau
                tau_all(i+1) = new_tau;
                tau_accept(i+1) = 1;
            else
                tau_all(i+1) = tau_all(i);
                tau_accept(i+1) = 0;
            end
        end
        
        curr_tau = tau_all(i+1);
        
        %%%%% Propose a new alpha %%%%%
        new_alpha = alpha_all(i) + alpha_epsilon * normrnd(0, 1);
        
        if new_alpha < min_alpha || new_alpha > max_alpha
            alpha_all(i+1) = alpha_all(i);
            alpha_accept(i+1) = 0;
        else
            log_alpha = compute_log_g(lambda, phi, curr_xi, curr_tau, new_alpha, gamma, label_data) ...
            - compute_log_g(lambda, phi, curr_xi, curr_tau, curr_alpha, gamma, label_data);
            transition_alpha = exp(log_alpha);
            if rand(1) < transition_alpha
                alpha_all(i+1) = new_alpha;
                alpha_accept(i+1) = 1;
            else
                alpha_all(i+1) = alpha_all(i);
                alpha_accept(i+1) = 0;
            end
        end      
    end
end

function g = compute_log_g(lambda, phi, xi, tau, alpha, gamma, label_data)
    g = -compute_phi(gamma, label_data, compute_T(xi, tau, alpha, lambda, phi))-0.5*norm(xi)^2;
end

function l = compute_phi(gamma, label_data, u)
    sum = norm((sign(u)-label_data).*abs(label_data))^2;
    l = sum/(2*gamma^2);
end

function T = compute_T(xi, tau, alpha, lambda, phi)
    x = (lambda + tau^2).^(-alpha/2) .* xi;
    T = convert_std_basis(x, phi);
end

function x = compute_rand_xi(num_data)
    x = normrnd(0, 1, num_data, 1);
end

function u = convert_std_basis(x, phi)
    u = phi*x;
end
