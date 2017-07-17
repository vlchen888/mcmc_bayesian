function [p, tau_mean, alpha_mean, M_mean] = test_mcmc_t_a_M(percent_fidelity, sigma)

    params = containers.Map;

    params('data_set') = string('mnist');
    params('laplacian') = string('self tuning');
    
    if params('data_set') == string('moons')
        N = 2000;
        data = moondata(1,100,N,sigma);
        params('data') = data;
        params('truth') = [-ones(floor(N/2)+1,1); ones(N-(floor(N/2)+1),1)];
        params('label_data') = generate_fidelity(percent_fidelity, params('truth'), length(data));
    elseif params('data_set') == string('voting')
        load('data3.mat')
        data = X;
        params('data') = data;
        params('truth') = [-ones(267,1); ones(168,1)];
        params('label_data') = generate_fidelity(percent_fidelity, params('truth'), length(data));
    elseif params('data_set') == string('mnist')
        load('mnist49data.mat')
        params('data') = data;
        load('mnist49truth.mat')
        truth = truth(:, 1) - truth(:, 2);
        params('truth') = truth;
        params('label_data') = generate_fidelity(percent_fidelity, truth, length(data));
    end
    
    figure(1)
    plot(params('label_data'))

    params('num_iterations') = 100000;
    burn_in = 5000;
    params('burn_in') = burn_in;
    params('movie') = true;

    % not used
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;
    %
    
    params('gamma') = 0.1;
    params('B') = 0.05;
    params('init_tau') = 1;
    params('init_alpha') = 35;
    params('init_M') = 5;
    
    params('min_tau')       = 0.1;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;
    
    params('min_M')       = 1;
    params('max_M')       = 300;
    
    params('alpha_epsilon') = 0;
    params('tau_epsilon')   = 0;
    
    params('max_M_jump') = 7;
    
        
    [tau_all, alpha_all, M_all, std, ~, ~, ~, ~] ...
        = mcmc_learn_t_a_M_noncentered(params);
    u_avg = mean(sign(std(:, burn_in:end)), 2); %avg the rows
        
    p = count_correct(u_avg, params('label_data'), params('truth'));
    
    tau_mean = mean(tau_all(burn_in:end));
    alpha_mean = mean(alpha_all(burn_in:end));
    M_mean = mean(M_all(burn_in:end));
end