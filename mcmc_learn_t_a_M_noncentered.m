function [tau_all, alpha_all, M_all, std, xi_accept, tau_accept, alpha_accept, M_accept] =...
        mcmc_learn_t_a_M_noncentered(params)
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
    init_M = params('init_M');
    
    min_tau = params('min_tau');
    max_tau = params('max_tau');
    min_alpha = params('min_alpha');
    max_alpha = params('max_alpha');
    min_M = params('min_M');
    max_M = params('max_M');
    
    alpha_epsilon = params('alpha_epsilon');
    tau_epsilon = params('tau_epsilon');
    
    if params('laplacian') == string('self tuning')
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == string('un')
        L = compute_laplacian_standard(data, p, q, l);
    end
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    [num_data, ~] = size(data);    
    xi_all = zeros(num_data, num_iterations);
    
    %%%%% Initialization from Fiedler Vector?? %%%%%
    xi_all(2, 1) = (lambda(2)+init_tau^2)^(init_alpha/2);
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
    
    M_all = zeros(1, num_iterations);
    M_all(1) = init_M;
    
    %%%%% Acceptance probabilities %%%%%
    xi_accept       = zeros(1, num_iterations);
    tau_accept      = zeros(1, num_iterations);
    alpha_accept    = zeros(1, num_iterations);
    M_accept        = zeros(1, num_iterations);
    
    %%%%% Store standard basis vectors %%%%%
    std = zeros(num_data, num_iterations);
        
    for i=1:num_iterations-1
        %%%%% Propose new state for xi %%%%%
        
        curr_xi = xi_all(:,i);
        curr_tau = tau_all(i);
        curr_alpha = alpha_all(i);
        curr_M = M_all(i);
        
        std(:, i) = compute_T(curr_xi, curr_tau, curr_alpha, curr_M, lambda, phi);
        
        x = compute_rand_xi(num_data);
        new_xi = (1-B^2)^0.5*curr_xi+B*x;
        u_curr = compute_T(curr_xi, curr_tau, curr_alpha, curr_M, lambda, phi);
        u_new = compute_T(new_xi, curr_tau, curr_alpha, curr_M, lambda, phi);
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
        if tau_epsilon > 0
            new_tau = curr_tau + tau_epsilon * normrnd(0, 1);
            if new_tau < min_tau || new_tau > max_tau
                tau_all(i+1) = tau_all(i);
                tau_accept(i+1) = 0;
            else
                log_tau = compute_log_g(lambda, phi, curr_xi, new_tau, curr_alpha, curr_M, gamma, label_data) ...
                - compute_log_g(lambda, phi, curr_xi, curr_tau, curr_alpha, curr_M, gamma, label_data);
                transition_tau = exp(log_tau);
                if rand(1) < transition_tau
                    tau_all(i+1) = new_tau;
                    tau_accept(i+1) = 1;
                else
                    tau_all(i+1) = tau_all(i);
                    tau_accept(i+1) = 0;
                end
            end
        else
            tau_all(i+1) = tau_all(i);
        end
        curr_tau = tau_all(i+1);
        
        %%%%% Propose a new alpha %%%%%
        if alpha_epsilon > 0
            new_alpha = curr_alpha + alpha_epsilon * normrnd(0, 1);

            if new_alpha < min_alpha || new_alpha > max_alpha
                alpha_all(i+1) = alpha_all(i);
                alpha_accept(i+1) = 0;
            else
                log_alpha = compute_log_g(lambda, phi, curr_xi, curr_tau, new_alpha, curr_M, gamma, label_data) ...
                - compute_log_g(lambda, phi, curr_xi, curr_tau, curr_alpha, curr_M, gamma, label_data);
                transition_alpha = exp(log_alpha);
                if rand(1) < transition_alpha
                    alpha_all(i+1) = new_alpha;
                    alpha_accept(i+1) = 1;
                else
                    alpha_all(i+1) = alpha_all(i);
                    alpha_accept(i+1) = 0;
                end
            end
        else
            alpha_all(i+1) = alpha_all(i);
        end
        curr_alpha = alpha_all(i+1);
        
        %%%%% Propose a new M %%%%%
        
        new_M = curr_M + compute_rndjump_M(params('max_M_jump'));
        
        if new_M < min_M || new_M > max_M
            M_all(i+1) = M_all(i);
            M_accept(i+1) = 0;
        else
            log_M = compute_log_g(lambda, phi, curr_xi, curr_tau, curr_alpha, new_M, gamma, label_data) ...
                - compute_log_g(lambda, phi, curr_xi, curr_tau, curr_alpha, curr_M, gamma, label_data);
            transition_M = exp(log_M);
            if rand(1) < transition_M
                M_all(i+1) = new_M;
                M_accept(i+1) = 1;
            else
                M_all(i+1) = M_all(i);
                M_accept(i+1) = 0;
            end
        end
        
    end
end

function g = compute_log_g(lambda, phi, xi, tau, alpha, M, gamma, label_data)
    g = -compute_phi(gamma, label_data, compute_T(xi, tau, alpha, M, lambda, phi))-0.5*norm(xi)^2;
end

function l = compute_phi(gamma, label_data, u)
    sum = norm((sign(u)-label_data).*abs(label_data))^2;
    l = sum/(2*gamma^2);
end

function T = compute_T(xi, tau, alpha, M, lambda, phi)
    mask = zeros(length(lambda), 1);
    mask(1:M) = 1;
    x = ((lambda + tau^2).^(-alpha/2) .* xi).*mask;
    T = convert_std_basis(x, phi);
end

function x = compute_rand_xi(num_data)
    x = normrnd(0, 1, num_data, 1);
end

function u = convert_std_basis(x, phi)
    u = phi*x;
end

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
