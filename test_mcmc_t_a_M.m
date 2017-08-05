function [p, M_mean, M_median, M_min] = test_mcmc_t_a_M(percent_fidelity, sigma)

    params = containers.Map;

    params('data_set') = "moons";
    params('laplacian') = "self tuning";
    
    if params('data_set') == "moons"
        N = 2000;
        data = moondata(1,100,N,sigma);
        params('data') = data;
        params('truth') = [-ones(floor(N/2)+1,1); ones(N-(floor(N/2)+1),1)];
        params('label_data') = generate_fidelity(percent_fidelity, params('truth'), length(data));
    elseif params('data_set') == "voting"
        load('data3.mat')
        data = X;
        params('data') = data;
        params('truth') = [-ones(267,1); ones(168,1)];
        params('label_data') = generate_fidelity(percent_fidelity, params('truth'), length(data));
    elseif params('data_set') == "mnist"
        load('mnist49data.mat')
        params('data') = data;
        load('mnist49truth.mat')
        truth = truth(:, 1) - truth(:, 2);
        params('truth') = truth;
        params('label_data') = generate_fidelity(percent_fidelity, truth, length(data));
    end
    
    params('num_iterations') = 100000;
    burn_in = 10000;
    params('burn_in') = burn_in;
    params('movie') = false;

    % not used
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;
    %
    
    params('gamma') = sigma;
    params('B') = 0.1;
    params('init_tau') = 2;
    params('init_alpha') = 35;
    params('init_M') = 50;
    
    params('min_tau')       = 0.1;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;
    
    params('min_M')       = 1;
    params('max_M')       = 70;
    
    params('alpha_epsilon') = 0;
    params('tau_epsilon')   = 0;
    
    params('max_M_jump') = 20;
    
        
    [M_all, sign_avg] = mcmc_learn_t_a_M_noncentered(params);
    p = count_correct(sign(sign_avg), params('label_data'), params('truth'));
    M_mean = mean(M_all(burn_in:end));
    M_median = median(M_all(burn_in:end));
    M_min = min(M_all(burn_in:end));
    
end