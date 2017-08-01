function [std] =...
        mcmc_learn_v_M(params)
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
    tau = params('tau');
    alpha = params('alpha');
    min_M = params('min_M');
    max_M = params('max_M');
    
    if params('laplacian') == "self tuning"
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == "un"
        L = compute_laplacian_standard(data, p, q, l);
    end
    
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    phi = phi(:, 1:max_M);
    lambda = lambda(1:max_M);
    
    % Arrays min_v, max_v
    min_v = (1-a) * (lambda + tau^2).^(-alpha/2);
    max_v = (1+a) * (lambda + tau^2).^(-alpha/2);
    
    std = zeros(length(data), num_iterations);
    v_all   = zeros(params('max_M'), num_iterations);
    v_all(:,1) = (lambda + tau^2).^(-alpha/2);
    v_accept = zeros(1, num_iterations);
    xi_all  = zeros(params('max_M'), num_iterations);
    xi_accept = zeros(1, num_iterations);
    M_all = zeros(1,num_iterations);
    M_all(1) = params('init_M');
    M_accept = zeros(1, num_iterations);
    
    tic;
    step_num = 1;
    for k=1:num_iterations-1
        curr_M = M_all(k);
        std(:, k) = compute_T(v_all(:,k), xi_all(:,k), curr_M, phi);
        
        if k >= params('burn_in') && mod(k,2500)==0
            curr_avg = mean(std(:,params('burn_in'):k), 2);
            fprintf('Sample number: %d, Elapsed time: %.4f\n', k, toc);
            fprintf('Classification accuracy: %f\n', count_correct(curr_avg, params('label_data'), params('truth')));
            fprintf('Acceptance rates: \n');
            fprintf('\txi: %.4f\n', mean(xi_accept(1:k)));
            fprintf('\tv: %.4f\n', mean(v_accept(1:k)));
            fprintf('\tM: %.4f\n', mean(M_accept(1:k)));
            
            figure(1)
            set(gcf, 'Position', [100, 300, 1000, 500])
            subplot(231)
            scatter_twomoons_classify(data, curr_avg, params('label_data'))
            title("Average u")
            
            subplot(232)
            plot(1:k, M_all(1:k), 1:k, movmean(M_all(1:k), [k 0]))
            legend('M trace', 'M running average');
            
            subplot(233)
            plot(movmean(xi_accept(1:k), [k 0]));
            title('\xi acceptance probability')
            
            subplot(234)
            plot(1:max_M, min_v, 1:max_M, max_v, 1:max_M, v_all(:,k))
            hold on
            line([curr_M curr_M],[0 max_v(1)], 'LineStyle', '--', 'Color', 'black');
            hold off
            legend('v_{min}', 'v_{max}', 'v_{curr}')
            
            subplot(235)
            plot(mean(v_all(:,params('burn_in'):k).*xi_all(:,params('burn_in'):k), 2));
            title("Average u_j")            
            
            subplot(236)
            histogram(M_all(1:k),'BinWidth',1);
            title('M histogram')
            
            drawnow
            
            fname = sprintf('figs/step_%i.png',step_num);
            print('-r144','-dpng',fname);
            step_num = step_num + 1;
            pause(.5);
        end
        
        % xi proposal
        hat_xi = (1-B^2)^0.5 * xi_all(:, k) + B*normrnd(0, 1, max_M, 1);
        log_xi = compute_phi(gamma, label_data, compute_T(v_all(:,k),xi_all(:,k),curr_M,phi)) - ...
            compute_phi(gamma, label_data, compute_T(v_all(:,k),hat_xi,curr_M,phi));
        if rand(1) < exp(log_xi)
            xi_all(:, k+1) = hat_xi;
            xi_accept(k+1) = 1;
        else
            xi_all(:, k+1) = xi_all(:, k);
            xi_accept(k+1) = 0;
        end
        
        % v proposal
        % could do non-identity diagonal matrix
        if epsilon > 0
            hat_v = v_all(:, k) + epsilon*normrnd(0,1,max_M,1);
            if ~isempty (find(hat_v > max_v, 1)) || ~isempty(find(hat_v < min_v, 1))
                v_all(:, k+1) = v_all(:, k);
                v_accept(k+1) = 0;
            else
                log_v = compute_phi(gamma, label_data, compute_T(v_all(:,k),xi_all(:,k+1),curr_M,phi)) - ...
                    compute_phi(gamma, label_data, compute_T(hat_v,xi_all(:,k+1),curr_M,phi));
                if rand(1) < exp(log_v)
                    v_all(:, k+1) = hat_v;
                    v_accept(k+1) = 1;
                else
                    v_all(:, k+1) = v_all(:, k);
                    v_accept(k+1) = 0;
                end
            end
        else
            v_all(:, k+1) = v_all(:, k);
            v_accept(k+1) = 1;
        end
        
        curr_v = v_all(:,k+1);
        curr_xi = xi_all(:,k+1);
        % M proposal
        M_hat = curr_M + compute_rndjump_M(params('max_M_jump'));
        if M_hat > max_M || M_hat < min_M
            M_all(k+1) = M_all(k);
            M_accept(k) = 0;
        else
            log_M = compute_phi(gamma, label_data, compute_T(curr_v, curr_xi, curr_M,phi)) - ...
                compute_phi(gamma, label_data, compute_T(curr_v, curr_xi, M_hat,phi));
            if rand(1) < exp(log_M)
                M_all(k+1) = M_hat;
                M_accept(k+1) = 1;
            else
                M_all(k+1) = M_all(k);
                M_accept(k) = 0;
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

function T = compute_T(v, xi, M, phi)
    T = phi(:,1:M)*(v(1:M).*xi(1:M));
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
