function [u_all, u_accept] = mcmc_gamma(params)
    num_iterations = params('num_iterations');
    
    data = params('data');
    label_data = params('label_data');
    p = params('p');
    q = params('q');
    l = params('l');
    tau = params('init_tau');
    alpha = params('init_alpha');
    gamma = params('gamma');
    B = params('B');
    
    if params('laplacian') == string('self tuning')
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == string('un')
        L = compute_laplacian_standard(data, p, q, l);
    end
    
    [num_data, ~] = size(data);

    lambda = eig(L);
    [phi, ~] = eig(L);
    M = 150;
    lambda = lambda(1:M);
    phi = phi(:,1:M);
    
            
    u_all = zeros(num_data, num_iterations);
    u_accept = zeros(1, num_iterations);
    
    %u(:, 1) = (lambda(1)+tau^2)^(-alpha/2)*phi(:,1);
    for i=1:num_iterations-1
        
        %%%% Propose new state for U %%%%
        u_star = (1-B^2)^0.5*u_all(:,i)+B*compute_rand(lambda, phi, tau, alpha);

        %%%%% COMPUTE TRANSITION PROBABILITY %%%%%
        log_uv = compute_log_likelihood(gamma, label_data, u_star) - ...
            compute_log_likelihood(gamma, label_data, u_all(:, i));
        transition_uv = exp(log_uv);

        if rand(1) < transition_uv %%%% TRANSITION %%%%
            u_all(:, i+1) = u_star;
            u_accept(i+1) = 1;
        else
            u_all(:, i+1) = u_all(:, i);
            u_accept(i+1) = 0;
        end
    end
end

function l = compute_log_likelihood(gamma, label_data, u)
    sum = norm((sign(u)-label_data).*abs(label_data))^2;
    l = -sum/(2*gamma^2);
end

function u = compute_rand(lambda, phi, tau, alpha)
    u = phi * ((eigvalweight(lambda, tau, alpha) .* normrnd(0, 1, length(lambda), 1)));
end

function w = eigvalweight(lambda, tau, alpha)
    w = (lambda + tau^2).^(-1/2*alpha);
end

