function [tau_all, alpha_all, std, u_accept, tau_accept, alpha_accept] =...
        mcmc_learn_t_a(Z, num_iterations, label_data, p, q, l, ...
        gamma, B, init_tau, init_alpha, min_tau, max_tau, min_alpha, ...
        max_alpha, alpha_epsilon, tau_epsilon)
    L = compute_laplacian_selftuning(Z);
    % [L, ~, ~] = compute_laplacian_standard(Z, p, q, l);
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    [num_data, ~] = size(Z);
        
    U = zeros(num_data, num_iterations);
    
    %%%%% Indicates initialization from Fiedler Vector %%%%%
    U(2, 1) = 1;
            
    tau_all = zeros(1, num_iterations);
    tau_all(1) = init_tau;
    
    alpha_all = zeros(1, num_iterations);
    alpha_all(1) = init_alpha;
    
    %%%%% Acceptance probabilities %%%%%
    u_accept        = zeros(1, num_iterations);
    tau_accept      = zeros(1, num_iterations);
    alpha_accept    = zeros(1, num_iterations);
    
    %%%%% Store standard basis vectors %%%%%
    std = zeros(num_data, num_iterations);
    std(:, 1) = convert_std_basis(U(:, 1), phi);
        
    for i=1:num_iterations-1
        
        %%%%% Propose new state for U %%%%%
        tau = tau_all(i);
        alpha = alpha_all(i);
        x = compute_rand_in_eigenbasis(lambda, tau, alpha);
        V_eigenbasis = (1-B^2)^0.5*U(:,i)+B*x;
        V = convert_std_basis(V_eigenbasis, phi);
        
        %%%%% Compute Transition Probability U -> V %%%%%
        log_uv = compute_log_likelihood(gamma, label_data, V) - ...
            compute_log_likelihood(gamma, label_data, std(:, i));
        transition_uv = exp(log_uv);
        
        %%%% Do transition %%%%
        if rand(1) < transition_uv 
            U(:, i+1) = V_eigenbasis;
            std(:, i+1) = V;
            u_accept(i+1) = 1;
        else
            U(:, i+1) = U(:, i);
            std(:, i+1) = std(:, i);
            u_accept(i+1) = 0;
        end
                
        %%%%% Propose a new tau %%%%%
        new_tau = tau_all(i) + tau_epsilon * normrnd(0, 1);
        if new_tau < min_tau
            new_tau = min_tau;
        elseif new_tau > max_tau
            new_tau = max_tau;
        end
        
        log_tau = compute_log_prior(U(:, i+1), lambda, new_tau, alpha_all(i)) ...
            - compute_log_prior(U(:, i+1), lambda, tau_all(i), alpha_all(i));
        transition_tau = exp(log_tau);
        if rand(1) < transition_tau
            tau_all(i+1) = new_tau;
            tau_accept(i+1) = 1;
        else
            tau_all(i+1) = tau_all(i);
            tau_accept(i+1) = 0;
        end
        
        %%%%% Propose a new alpha %%%%%
        new_alpha = alpha_all(i) + alpha_epsilon * normrnd(0, 1);
        
        if new_alpha < min_alpha
            new_alpha = min_alpha;
        elseif new_alpha > max_alpha
            new_alpha = max_alpha;
        end
                
        log_alpha = compute_log_prior(U(:, i+1), lambda, tau_all(i+1), new_alpha) ...
            - compute_log_prior(U(:, i+1), lambda, tau_all(i+1), alpha_all(i));
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

function log_prior = compute_log_prior(u, lambda, tau, alpha)
    log_prior = -0.5*compute_log_det_c(lambda, tau, alpha) ...
        -0.5*compute_inner_prod(u, lambda, tau, alpha);
end

function innerProd = compute_inner_prod(u, lambda, tau, alpha)
    innerProd = 0;
    for j = 1:length(lambda)
        innerProd = innerProd + u(j)^2 * (lambda(j) + tau^2)^alpha;
    end
end

function log_det_c = compute_log_det_c(lambda, tau, alpha)
    sum = 0;
    for j = 1:length(lambda)
        sum = sum + log(lambda(j) + tau^2);
    end
    log_det_c = -alpha * sum;
end

function l = compute_log_likelihood(gamma, label_data, u)
    sum = 0;
    for i=1:length(label_data)
        %%%% Check if i \in Z' %%%%
        if label_data(i) ~= 0
            sum = sum + abs(sign(u(i)) - label_data(i))^2;
        end
    end
    l = -sum/(2*gamma^2);
end

function x = compute_rand_in_eigenbasis(lambda, tau, alpha)
    x = zeros(length(lambda), 1);
    for j=1:length(lambda)
       zeta = normrnd(0, 1);
       x(j) = eigvalweight(lambda(j), tau, alpha) * zeta;
    end
end

function u = convert_std_basis(x, phi)
    u = 0;
    for j=1:length(phi)
        u = u + x(j)* phi(:,j);
    end
end

function w = eigvalweight(lambda, tau, alpha)
    w = (lambda + tau^2)^(-alpha/2);
end