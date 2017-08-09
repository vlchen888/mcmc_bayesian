function u_avg = mcmc_multiclass(params)
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
    u_avg = zeros(num_data, k);
    sign_avg = zeros(num_data, k);
    
    %%%%% Acceptance probabilities %%%%%
    xi_accept = zeros(k, num_iterations);
    
    tic;    
    for i=1:num_iterations-1
        u_curr = phi * ((lambda + tau^2).^(-alpha/2) .* xi_curr);
        sign_curr = compute_S_multiclass(u_curr, k);
        if i >= params('burn_in')
            u_avg = (u_avg * (i-params('burn_in')) + u_curr)/(i-params('burn_in') + 1);
            sign_avg = (sign_avg * (i-params('burn_in')) + sign_curr)/(i-params('burn_in') + 1);
        end
        
        if i >= params('burn_in') && mod(i, 2500) == 0            
            p = count_correct_multiclass(compute_S_multiclass(u_avg, k), params('label_data'), params('truth'));
            q = count_correct_multiclass(compute_S_multiclass(sign_avg, k), params('label_data'), params('truth'));
            fprintf('Sample number: %d, Time elapsed: %.2f\n', i, toc);
            fprintf('Classification accuracy with S(E(u)): %.4f\n', p);
            fprintf('Classification accuracy with S(E(S(u))): %.4f\n', q);
            fprintf('\txi accept acceptance probability: %.4f\n', mean(xi_accept(:,1:i),2));
            figure(1)
            set(gcf, 'Position', [0, 500, 500, 400])
            mnist_heatmap(compute_S_multiclass(u_avg, k), params('truth'), params('digs'), "S(E(u)) classification");
            
            figure(2)
            set(gcf, 'Position', [0, 500, 500, 400])
            mnist_heatmap(compute_S_multiclass(sign_avg, k), params('truth'), params('digs'), "S(E(S(u))) classification");
            
            figure(3)
            set(gcf, 'Position', [0, 0, 600, 400])
            subplot(311)
            plot(u_curr)
            title('Current u')
            subplot(312)
            plot((lambda + tau^2).^(-alpha/2) .* xi_curr)
            title('Current u_j')
            subplot(313)
            plot(u_avg);
            title('Average u')
            
            figure(4)
            set(gcf, 'Position', [0, 0, 600, 400])
            for kk = 1:k
                subplot(k,1,kk)
                plot(movmean(xi_accept(kk,1:i),[i 0]))
                if kk == 1
                    title('\xi acceptance probability')
                end
            end
            
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
    diff = abs(compute_S_multiclass(u, k) - label_data)/sqrt(2);
    diff = diff(sum(label_data, 2) ~= 0 , :);
    l = sum(sum(diff))/(2*gamma^2);
end

function T = compute_T(xi, tau, alpha, lambda, phi)
    T = phi* ( (lambda + tau^2).^(-alpha/2) .* xi);
end