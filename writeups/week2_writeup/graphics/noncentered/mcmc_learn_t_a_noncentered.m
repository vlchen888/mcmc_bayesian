function [tau_all, alpha_all, std, xi_accept, tau_accept, alpha_accept] =...
        mcmc_learn_t_a_noncentered(Z, num_iterations, label_data, p, q, l, ...
        gamma, B, init_tau, init_alpha, min_tau, max_tau, min_alpha, ...
        max_alpha, alpha_epsilon, tau_epsilon)
    L = compute_laplacian_selftuning(Z);
    %[L, ~, ~] = computeLaplacian(Z, p, q, l);
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    [num_data, ~] = size(Z);
    
    xi_all = zeros(num_data, num_iterations);
    
    %%%%% Initialization from Fiedler Vector?? %%%%%
    %xi_all(2, 1) = (lambda(2)+init_tau^2)^(init_alpha/2);
    %xi_all(2, 1) = 1;
    
    %%%%% Noisy initialization?? %%%%%
    %xi_all(:,1) = compute_rand_xi(num_data)/1000;
    %xi_all(2, 1) = 1;
     
    %%%%% Random initialization? %%%%%
    %xi_all(:,1) = compute_rand_xi(num_data)/2;
    
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
        
    for i=1:num_iterations-1
        %%%%% Propose new state for xi %%%%%
        curr_xi = xi_all(:,i);
        curr_tau = tau_all(i);
        curr_alpha = alpha_all(i);
        
        std(:, i) = compute_T(curr_xi, curr_tau, curr_alpha, lambda, phi);
        
        x = compute_rand_xi(num_data);
        new_xi = (1-B^2)^0.5*curr_xi+B*x;
        
        log_xi_trans = compute_log_g(lambda, phi, new_xi, curr_tau, curr_alpha, gamma, label_data) ...
            - compute_log_g(lambda, phi, curr_xi, curr_tau, curr_alpha, gamma, label_data);
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
        if new_tau < min_tau
            new_tau = min_tau;
        elseif new_tau > max_tau
            new_tau = max_tau;
        end
        
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
        
        curr_tau = tau_all(i+1);
        
        %%%%% Propose a new alpha %%%%%
        new_alpha = alpha_all(i) + alpha_epsilon * normrnd(0, 1);
        
        if new_alpha < min_alpha
            new_alpha = min_alpha;
        elseif new_alpha > max_alpha
            new_alpha = max_alpha;
        end
                
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

function g = compute_log_g(lambda, phi, xi, tau, alpha, gamma, label_data)
    g = -compute_phi(gamma, label_data, compute_T(xi, tau, alpha, lambda, phi))-0.5*norm(xi)^2;
end

function l = compute_phi(gamma, label_data, u)
    sum = 0;
    for i=1:length(label_data)
        %%%% Check if i \in Z' %%%%
        if label_data(i) ~= 0
            sum = sum + abs(sign(u(i)) - label_data(i))^2;
        end
    end
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
    u = 0;
    for j=1:length(phi)
        u = u + x(j)* phi(:,j);
    end
end
