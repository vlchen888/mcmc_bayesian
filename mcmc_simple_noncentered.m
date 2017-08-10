function cont = mcmc_simple_noncentered(params)
    data = params('data');
    num_iterations = params('num_iterations');
    label_data = params('label_data');
    p = params('p');
    q = params('q');
    l = params('l');
    
    gamma = params('gamma');
    B = params('B');
    
    curr_tau = params('init_tau');
    curr_alpha = params('init_alpha');
    curr_M = params('init_M');
    
    %tic;
    if params('laplacian') == "self tuning"
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == "un"
        L = compute_laplacian_standard(data, p, q, l);
    end
    %toc
    
    %fprintf('Laplacian computation complete!\n');
    
    lambda = eig(L);
    [phi, ~] = eig(L);
    lambda = lambda(1:curr_M);
    phi = phi(:, 1:curr_M);
    
    [num_data, ~] = size(data);
    curr_xi = zeros(curr_M, 1);
    %%%%% Initialization from Fiedler Vector?? %%%%%
    %xi_all(2, 1) = (lambda(2)+init_tau^2)^(-init_alpha/2);
    
    uj_avg = zeros(1, curr_M);
        
    %%%%% Acceptance probabilities %%%%%
    xi_accept       = zeros(1, num_iterations);
    
    %%%%% Store standard basis vectors %%%%%
    sign_avg = zeros(num_data, 1);
    
    xi_all = zeros(curr_M, num_iterations);
    v_all = zeros(curr_M, num_iterations);
    
    %%%%% Store correct percent %%%%%
    correct_p = zeros(1, num_iterations);
    tic;
    for i=1:num_iterations-1        
        xi_all(:,i) = curr_xi;
        v_all(:,i) = [0; (lambda(2:end) + curr_tau^2).^(-curr_alpha/2)/norm_const(lambda,curr_tau,curr_alpha)];
        curr_u = compute_T(curr_xi, curr_tau, curr_alpha, curr_M, lambda, phi);
        if i>= params('burn_in')
            sign_avg = (sign(curr_u) + sign_avg*(i-params('burn_in')))/(i-params('burn_in')+1);
            uj_avg = ( v_all(:,i).*xi_all(:,i) + uj_avg*(i-params('burn_in')))/(i-params('burn_in')+1);
        end
        
        %%%% Movie things %%%%
        if params('movie') && i >= params('burn_in') && mod(i,2500)==0
            curr_avg = sign_avg;
            p = count_correct(curr_avg, params('label_data'), params('truth'));
            correct_p(i) = p;
            
            fprintf('Sample number: %d, Elapsed time: %.4f\n', i, toc);
            fprintf('Classification accuracy S(E(S(u))): %f\n', p);
            figure(2)
            subplot(231)
            set(gcf, 'Position', [100, 300, 1000, 500])
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


            subplot(232)
            plot(movmean(xi_accept(1:i), [i 0]));
            title('\xi acceptance probability')
            
            subplot(233)
            plot(uj_avg);
            title('Average u_j')
            
            subplot(236)
            plot(correct_p(correct_p~=0))
            title('Classification accuracy');

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
    end
    
    cont = containers.Map;
    cont('v_all') = v_all;
    cont('xi_all') = xi_all;
end

function l = compute_phi(gamma, label_data, u)
    sum = norm((sign(u)-label_data).*abs(label_data))^2;
    l = sum/(2*gamma^2);
end

function T = compute_T(xi, tau, alpha, M, lambda, phi)
    x = (lambda(2:M) + tau^2).^(-alpha/2) .* xi(2:M);
    T = phi(:,2:M)*x;
end

function x = compute_rand_xi(num_data)
    x = normrnd(0, 1, num_data, 1);
end

function const = norm_const(lambda, tau, alpha)
    const = sqrt(sum((lambda(2:end)+tau^2).^-alpha))/sqrt(length(lambda(2:end)));
end
