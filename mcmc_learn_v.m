function [std, v_all, xi_all, v_accept, xi_accept] =...
        mcmc_learn_v(params)
    data = params('data');
    num_iterations = params('num_iterations');
    label_data = params('label_data');
    p = params('p');
    q = params('q');
    l = params('l');
    
    gamma = params('gamma');
    B = params('B');
    a = params('a');
    epsilon = params('epsilon');
    
    if params('laplacian') == string('self tuning')
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == string('un')
        L = compute_laplacian_standard(data, p, q, l);
    end
    
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    M = 50;
    phi = phi(:, 1:M);
    lambda = lambda(1:M);
    
    % Arrays min_v, max_v
    load('moon_E_u_sq.mat');
    min_v = (1-a) * E_u_sq.^0.5;
    max_v = (1+a) * E_u_sq.^0.5;
    
    std = zeros(length(data), num_iterations);
    v_all   = zeros(M, num_iterations);
    v_all(:,1) = E_u_sq.^0.5;
    v_accept = zeros(1, num_iterations);
    xi_all  = zeros(M, num_iterations);
    xi_accept = zeros(1, num_iterations);
    
    for k=1:num_iterations-1
        std(:, k) = compute_T(v_all(:,k), xi_all(:,k), phi);
        hat_xi = (1-B^2)^0.5 * xi_all(:, k) + B*normrnd(0, 1, M, 1);
        log_xi = compute_phi(gamma, label_data, compute_T(v_all(:,k),xi_all(:,k),phi)) - ...
            compute_phi(gamma, label_data, compute_T(v_all(:,k),hat_xi,phi));
        if rand(1) < exp(log_xi)
            xi_all(:, k+1) = hat_xi;
            xi_accept(k+1) = 1;
        else
            xi_all(:, k+1) = xi_all(:, k);
            xi_accept(k+1) = 0;
        end
        
        % could do non-identity diagonal matrix
        hat_v = v_all(:, k) + epsilon*normrnd(0,1,M,1);
        if ~isempty (find(hat_v > max_v, 1)) || ~isempty(find(hat_v < min_v, 1))
            v_all(:, k+1) = v_all(:, k);
            v_accept(k+1) = 0;
        else
            log_v = compute_phi(gamma, label_data, compute_T(v_all(:,k),xi_all(:,k+1),phi)) - ...
                compute_phi(gamma, label_data, compute_T(hat_v,xi_all(:,k+1),phi));
            if rand(1) < exp(log_v)
                v_all(:, k+1) = hat_v;
                v_accept(k+1) = 1;
            else
                v_all(:, k+1) = v_all(:, k);
                v_accept(k+1) = 0;
            end
        end
    end
end

function h = compute_log_h(v, xi, phi, label_data, gamma)
    h = -compute_phi(gamma, label_data, T(v, xi, phi)) - 0.5*norm(xi)^2;
end

function l = compute_phi(gamma, label_data, u)
    sum = norm((sign(u)-label_data).*abs(label_data))^2;
    l = sum/(2*gamma^2);
end

function T = compute_T(v, xi, phi)
    T = phi*(v.*xi);
end
