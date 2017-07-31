function [signs_all] = mcmc_multiclass(params)
    data = params('data');
    num_iterations = params('num_iterations');
    label_data = params('label_data');
    p = params('p');
    q = params('q');
    l = params('l');
    gamma = params('gamma');
    B = params('B');
    tau = params('tau');
    alpha = params('alpha');
    k = params('k');
    
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
    
    xi_curr = zeros(M, k);
    signs_all = zeros(k, num_data, num_iterations);
    
    %%%%% Acceptance probabilities %%%%%
    xi_accept       = zeros(k, num_iterations);
    
    tic;    
    for i=1:num_iterations-1
        signs_all(:,:,i) = compute_S(compute_T(xi_curr(:,:), tau, alpha, lambda, phi), k);
        
        if mod(i, 2500) == 0
            figure(1)
            avg_label = mean(signs_all(:,:,1:i), 3);
            [~, I] = max(avg_label);
            A = eye(k);
            curr_label = A(:,I);
            p = count_correct_multiclass(curr_label, params('label_data'), params('truth'));
            fprintf('Sample number: %d, Time elapsed: %.2f\n', i, toc);
            fprintf('Classification accuracy: %.4f\n', p);
            fprintf('\txi accept acceptance probability: %.4f\n', mean(xi_accept(:,1:i),2));
            mnist_heatmap(curr_label, params('truth'), params('digs'));
            drawnow
        end
        
        for j=1:k
            xi_hat = sqrt(1-B^2)*xi_curr(:,j) + B*normrnd(0,1,M,1);
            xi_new = xi_curr;
            xi_new(:, j) = xi_hat;
            log_xi = compute_phi(gamma, label_data, compute_T(xi_curr, tau, alpha, lambda, phi), k) - ...
                compute_phi(gamma, label_data, compute_T(xi_new, tau, alpha, lambda, phi), k);
            if rand(1) < exp(log_xi)
                xi_curr(:,j) = xi_hat;
                xi_accept(j,i+1) = 1;
            else
                xi_accept(j,i+1) = 0;
            end
        end
    end
    
end

function l = compute_phi(gamma, label_data, u, k)
    diff = abs(compute_S(u, k) - label_data)/sqrt(2);
    diff = diff(:, find (sum(label_data) ~= 0) );
    l = sum(sum(diff))/(2*gamma^2);
end

function S = compute_S(u, k)
    [~, I] = max(u,[],2);
    A = eye(k);
    S = A(:, I);
end

function T = compute_T(xi, tau, alpha, lambda, phi)
    T = phi* ( (lambda + tau^2).^(-alpha/2) .* xi);
end