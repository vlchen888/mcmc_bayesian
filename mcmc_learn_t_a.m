function [tau_all, alpha_all, std, u_accept, tau_accept, alpha_accept] =...
        mcmc_learn_t_a(params)
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
    
    M = 50;
    
    if params('laplacian') == string('self tuning')
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == string('un')
        L = compute_laplacian_standard(data, p, q, l);
    end
    lambda = eig(L);
    [phi, ~] = eig(L);
    phi = phi(:, 1:M);
    lambda = lambda(1:M);
    
    [num_data, ~] = size(data);
    
    U = zeros(M, num_iterations);
    
    %%%%% Indicates initialization from Fiedler Vector %%%%%
    U(1, 1) = (lambda(1) + init_tau^2)^(-init_alpha/2);
    U(2, 1) = -(lambda(2) + init_tau^2)^(-init_alpha/2);
            
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
        if new_tau < min_tau || new_tau > max_tau
            tau_all(i+1) = tau_all(i);
            tau_accept(i+1) = 0;
        else
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
        end
        
        %%%%% Propose a new alpha %%%%%
        new_alpha = alpha_all(i) + alpha_epsilon * normrnd(0, 1);
        
        if new_alpha < min_alpha || new_alpha > max_alpha
            alpha_all(i+1) = alpha_all(i);
            alpha_accept(i+1) = 0;
        else
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
end

function log_prior = compute_log_prior(u, lambda, tau, alpha)
    log_prior = -0.5*compute_log_det_c(lambda, tau, alpha) ...
        -0.5*compute_inner_prod(u, lambda, tau, alpha);
end

function innerProd = compute_inner_prod(u, lambda, tau, alpha)
    innerProd = (u.^2)' * (lambda + tau^2).^alpha;
end

function log_det_c = compute_log_det_c(lambda, tau, alpha)
    log_det_c = -alpha * sum(log(lambda + tau^2));
end

function l = compute_log_likelihood(gamma, label_data, u)
    sum = norm((sign(u)-label_data).*abs(label_data))^2;
    l = -sum/(2*gamma^2);
end

function x = compute_rand_in_eigenbasis(lambda, tau, alpha)
    x = eigvalweight(lambda, tau, alpha) .* normrnd(0, 1, length(lambda), 1);
end

function u = convert_std_basis(x, phi)
    u = phi*x;
end

function w = eigvalweight(lambda, tau, alpha)
    w = (lambda + tau^2).^(-alpha/2);
end