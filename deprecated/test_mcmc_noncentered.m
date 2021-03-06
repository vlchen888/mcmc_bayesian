function [p, tau_mean, alpha_mean] = test_mcmc_hier_noncentered(percent_fidelity, sigma)
    params = containers.Map;

    params('data_set') = "mnist";
    params('laplacian') = "self tuning";
    
    if params('data_set') == "moons"
        N = 2000;
        data = moondata(1,100,N,sigma);
        params('data') = data;
        params('label_data') = generate_moons_fidelity(percent_fidelity, N);
        params('truth') = [zeros(floor(N/2)+1,1) - 1; zeros(N-(floor(N/2)+1),1) + 1];
    elseif params('data_set') == "mnist"
        load('mnist49data.mat')
        params('data') = data;
        load('mnist49truth.mat')
        truth = truth(:,1) - truth(:,2);
        params('truth') = truth;
        params('label_data') = generate_fidelity(percent_fidelity, params('truth'), length(data));
    end

    params('num_iterations') = 100000;
    burn_in = 5000;

    params('p') = 2;
    params('q') = 2;
    params('l') = 1;

    params('gamma')         = 0.1;
    params('B')             = 0.1;
    params('init_tau')      = 2;
    params('init_alpha')    = 1;

    params('min_tau')       = 0.01;
    params('max_tau')       = 35;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;

    params('alpha_epsilon') = 0;
    params('tau_epsilon')   = 0;


    [tau_all, alpha_all, std, xi_accept, tau_accept, alpha_accept] = mcmc_learn_t_a_noncentered(params);

    u_avg = mean(std(:, burn_in:end), 2); %avg the rows
    
    if params('data_set') == "moons"
        figure(1)
        clf
        scatter_twomoons_classify(data, u_avg, params('label_data'))

        tau_mean = mean(tau_all(burn_in:end));
        alpha_mean = mean(alpha_all(burn_in:end));
    elseif params('data_set') == "mnist"
        
    end
    
    p = count_correct(u_avg, params('label_data'), params('truth'));

    figure(2)
    clf
    subplot(2,1,1)
    plot(tau_all)
    xlabel('\tau')
    subplot(2,1,2)
    plot(alpha_all)
    xlabel('\alpha')

    xi_avg_accept = mean(xi_accept(burn_in:end))
    tau_avg_accept = mean(tau_accept(burn_in:end))
    alpha_avg_accept = mean(alpha_accept(burn_in:end))
end
