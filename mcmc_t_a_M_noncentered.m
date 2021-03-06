function cont = mcmc_t_a_M_noncentered(params)
    data = params('data');
    num_iterations = params('num_iterations');
    label_data = params('label_data');
    
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
    
    if params('laplacian') == "self tuning"
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == "un"
        p = params('p');
        q = params('q');
        l = params('l');
        L = compute_laplacian_standard(data, p, q, l);
    end
    
    lambda = eig(L);
    [phi, ~] = eig(L);
    lambda = lambda(1:params('max_M'));
    phi = phi(:, 1:params('max_M'));
    
    if params('remove-zero-eig')
        lambda(1) = Inf;
    end
    
    [num_data, ~] = size(data);
    curr_xi = zeros(params('max_M'), 1);
    
    tau_all = zeros(1, num_iterations);
    tau_all(1) = init_tau;
    
    alpha_all = zeros(1, num_iterations);
    alpha_all(1) = init_alpha;
    
    uj_avg = zeros(1, params('max_M'));
    
    M_all = zeros(1, num_iterations);
    M_all(1) = init_M;
    
    %%%%% Acceptance probabilities %%%%%
    xi_accept       = zeros(1, num_iterations);
    tau_accept      = zeros(1, num_iterations);
    alpha_accept    = zeros(1, num_iterations);
    M_accept        = zeros(1, num_iterations);
    
    %%%%% Store standard basis vectors %%%%%
    sign_avg = zeros(num_data, 1);
    
    xi_all = zeros(params('max_M'), num_iterations);
    v_all = zeros(params('max_M'), num_iterations);
    
    %%%%% Store correct percent %%%%%
    accuracy_map = containers.Map('KeyType','int32','ValueType','double')
    tic;
    for i=1:num_iterations-1
        curr_tau = tau_all(i);
        curr_alpha = alpha_all(i);
        curr_M = M_all(i);
        
        xi_all(:,i) = curr_xi;
        v_tmp = (lambda + curr_tau^2).^(-curr_alpha/2)/norm_const(lambda,curr_tau,curr_alpha);
        v_tmp(curr_M+1:end) = 0;
        v_all(:,i) = v_tmp;
        
        
        curr_u = compute_T(curr_xi, curr_tau, curr_alpha, curr_M, lambda, phi);
        if i>= params('burn_in')
            sign_avg = (sign(curr_u) + sign_avg*(i-params('burn_in')))/(i-params('burn_in')+1);
            uj_avg = ( v_all(:,i).*xi_all(:,i) + uj_avg*(i-params('burn_in')))/(i-params('burn_in')+1);
        end
        
        if mod(i, params('accuracy_often')) == 0
            accuracy_map(i) = count_correct(sign_avg, params('label_data'), params('truth'));
        end
        
        %%%% Movie things %%%%
        if params('movie') && i >= params('burn_in') && mod(i,params('movie_often'))==0
            curr_avg = sign_avg;
            p = count_correct(curr_avg, params('label_data'), params('truth'));
            
            fprintf('Sample number: %d, Elapsed time: %.4f\n', i, toc);
            fprintf('Classification accuracy S(E(S(u))): %f\n', p);
            fprintf('Acceptance rates: \n\txi:%f \n\ttau:%f \n\talpha:%f \n\tM:%f \n',...
            sum(xi_accept)/i,sum(tau_accept)/i,sum(alpha_accept)/i,sum(M_accept)/i);
            figure(2)
            subplot(221)
            set(gcf, 'Position', [0, 1000, 600, 500])
            if params('data_set') == "moons"
                scatter_twomoons_classify(data, curr_avg, params('label_data'))
                title('Average u scatter')
            elseif params('data_set') == "voting"
                plotBar(curr_avg);
                title('Average u')
            elseif params('data_set') == "mnist"
                plotBar(curr_avg);
                title('Average u')
            end

            subplot(222)
            plot(movmean(xi_accept(1:i), [i 0]));
            title('\xi acceptance probability')
            
            subplot(223)
            plot(uj_avg);
            title('Average u_j')
            
            subplot(224)
            plot(cell2mat(accuracy_map.keys), cell2mat(accuracy_map.values))
            title('Classification accuracy');
            xlabel('Sample number');
            
            figure(3)
            hyp_count = 0;
            if tau_epsilon > 0
                hyp_count = hyp_count + 1;
            end
            
            if alpha_epsilon > 0
                hyp_count = hyp_count + 1;
            end
            
            if params('max_M_jump') > 0
                hyp_count = hyp_count + 1;
            end
            set(gcf, 'Position', [700, 1000, 700, 200*hyp_count])
            fig_count = 1;
            
            if tau_epsilon > 0
                subplot(hyp_count,2,fig_count)
                plot(1:i,tau_all(1:i),1:i,movmean(tau_all(1:i), [i 0]));
                legend('\tau trace', '\tau running average');

                subplot(hyp_count,2,fig_count+1)
                plot(movmean(tau_accept(1:i),[i 0]))
                title('\tau acceptance probability')
                
                fig_count = fig_count + 2;
            end
            
            if alpha_epsilon > 0
                subplot(hyp_count,2,fig_count)
                plot(1:i,alpha_all(1:i),1:i,movmean(alpha_all(1:i), [i 0]));
                legend('\alpha trace', '\alpha running average');

                subplot(hyp_count,2,fig_count+1)
                plot(movmean(alpha_accept(1:i),[i 0]))
                title('\alpha acceptance probability')
                
                fig_count = fig_count + 2;
            end
            
            if params('max_M_jump') > 0
                subplot(hyp_count,2,fig_count)
                plot(1:i, M_all(1:i), 1:i, movmean(M_all(1:i), [i 0]));
                legend('M trace', 'M running average');

                subplot(hyp_count,2,fig_count+1)
                plot(movmean(M_accept(1:i), [i 0]));
                title('M acceptance probability')
                
                %fig_count = fig_count + 2;
            end

            drawnow
            pause(.5)

            %fname = sprintf('figs/step_%i.png',step_num);
            %print('-r144','-dpng',fname);
            %step_num = step_num + 1;
        end
        
        %%%%% Propose new state for xi %%%%%
        new_xi = (1-B^2)^0.5*curr_xi+B*compute_rand_xi(length(lambda));
        u_old = compute_T(curr_xi, curr_tau, curr_alpha, curr_M, lambda, phi);
        u_hat = compute_T(new_xi, curr_tau, curr_alpha, curr_M, lambda, phi);
        log_xi_trans = compute_phi(gamma, label_data, u_old) - compute_phi(gamma, label_data, u_hat);
        transition_xi = exp(log_xi_trans);
        
        %%%% Do transition %%%%
        if rand(1) < transition_xi
            curr_xi = new_xi;
            xi_accept(i+1) = 1;
        else
            xi_accept(i+1) = 0;
        end
                        
        %%%%% Propose a new tau %%%%%
        if tau_epsilon > 0
            tau_hat = curr_tau + tau_epsilon * normrnd(0, 1);
            if tau_hat < min_tau || tau_hat > max_tau
                tau_all(i+1) = tau_all(i);
                tau_accept(i+1) = 0;
            else
                u_old = compute_T(curr_xi, curr_tau, curr_alpha, curr_M, lambda, phi);
                u_hat = compute_T(curr_xi, tau_hat, curr_alpha, curr_M, lambda, phi);
                log_tau = compute_phi(gamma, label_data, u_old) - compute_phi(gamma, label_data, u_hat);
                transition_tau = exp(log_tau);
                if rand(1) < transition_tau
                    tau_all(i+1) = tau_hat;
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
            alpha_hat = curr_alpha + alpha_epsilon * normrnd(0, 1);

            if alpha_hat < min_alpha || alpha_hat > max_alpha
                alpha_all(i+1) = alpha_all(i);
                alpha_accept(i+1) = 0;
            else
                u_old = compute_T(curr_xi, curr_tau, curr_alpha, curr_M, lambda, phi);
                u_hat = compute_T(curr_xi, curr_tau, alpha_hat, curr_M, lambda, phi);
                log_alpha = compute_phi(gamma, label_data, u_old) - compute_phi(gamma, label_data, u_hat);
                transition_alpha = exp(log_alpha);
                if rand(1) < transition_alpha
                    alpha_all(i+1) = alpha_hat;
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
        if params('max_M_jump') > 0
            M_hat = curr_M + compute_rndjump_M(params('max_M_jump'));

            if M_hat < min_M || M_hat > max_M
                M_all(i+1) = M_all(i);
                M_accept(i+1) = 0;
            else
                u_old = compute_T(curr_xi, curr_tau, curr_alpha, curr_M, lambda, phi);
                u_hat = compute_T(curr_xi, curr_tau, curr_alpha, M_hat, lambda, phi);
                log_M = compute_phi(gamma, label_data, u_old) - compute_phi(gamma, label_data, u_hat);            
                transition_M = exp(log_M);
                if rand(1) < transition_M
                    M_all(i+1) = M_hat;
                    M_accept(i+1) = 1;
                else
                    M_all(i+1) = M_all(i);
                    M_accept(i+1) = 0;
                end
            end
        else
            M_all(i+1) = M_all(i);
        end
    end
    
    cont = containers.Map;
    cont('v_all') = v_all;
    cont('xi_all') = xi_all;
    cont('tau_all') = tau_all;
    cont('alpha_all') = alpha_all;
    cont('M_all') = M_all;
end

function l = compute_phi(gamma, label_data, u)
    sum = norm((sign(u)-label_data).*abs(label_data))^2;
    l = sum/(2*gamma^2);
end

function T = compute_T(xi, tau, alpha, M, lambda, phi)
    x = (lambda(1:M) + tau^2).^(-alpha/2) .* xi(1:M);
    T = phi(:,1:M)*x;
end

function x = compute_rand_xi(num_data)
    x = normrnd(0, 1, num_data, 1);
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

function const = norm_const(lambda, tau, alpha)
    const = sqrt(sum((lambda+tau^2).^-alpha))/sqrt(length(lambda));
end
