function p = test_mcmc_t_a_M(percent_fidelity, sigma)

    params = containers.Map;

    N = 2000;
    data = moondata(1,100,N,sigma);
    params('data') = data;
    params('label_data') = generate_moons_fidelity(percent_fidelity, N);

    params('data_set') = string('moons');
    params('laplacian') = string('self tuning');

    params('num_iterations') = 1000;
    burn_in = 100;

    % not used
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;
    %
    
    params('gamma') = 0.1;
    params('B') = 0.1;
    params('init_tau') = 1;
    params('init_alpha') = 1;
    params('init_M') = 70;
    
    params('min_tau')       = 0.1;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;
    
    params('min_M')       = 1;
    params('max_M')       = 1000;
    
    params('alpha_epsilon') = 0.1;
    params('tau_epsilon')   = 0.1;


    [tau_all, alpha_all, M_all, std, u_accept, tau_accept, alpha_accept, M_accept] ...
        = mcmc_learn_t_a_M(params);
u_avg = mean(sign(std(:, burn_in:end)), 2); %avg the rows

    figure(1)
    clf
    scatter_twomoons_classify(data, u_avg, params('label_data'))
    p = count_correct(u_avg, params('label_data'), [zeros(floor(N/2)+1,1) - 1; zeros(N-(floor(N/2)+1),1) + 1]);
    %tau_mean = mean(tau_all(burn_in:end));
    %alpha_mean = mean(alpha_all(burn_in:end));
    figure(2)
    clf
    subplot(3,1,1)
    plot(tau_all)
    xlabel('\tau');
    subplot(3,1,2)
    plot(alpha_all)
    xlabel('\alpha');
    subplot(3,1,3)
    plot(M_all)
    xlabel('M');
    
    u_avg_accept = mean(u_accept(burn_in:end))
    tau_avg_accept = mean(tau_accept(burn_in:end))
    alpha_avg_accept = mean(alpha_accept(burn_in:end))
    M_avg_accept = mean(M_accept(burn_in:end))

end