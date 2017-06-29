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
    
    if params('laplacian') == string('self tuning')
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == string('un')
        [L, ~, ~] = compute_laplacian_standard(data, p, q, l);
    end
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    [num_data, ~] = size(data);
    
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
