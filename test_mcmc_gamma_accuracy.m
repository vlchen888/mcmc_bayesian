function p = test_mcmc_gamma_accuracy(percent_fidelity, sigma)

    params = containers.Map;

    N = 2000;
    data = moondata(1,100,N,sigma);
    params('data') = data;
    params('label_data') = generate_moons_fidelity(percent_fidelity, N);

    params('data_set') = string('moons');
    params('laplacian') = string('self tuning');

    params('num_iterations') = 100000;
    burn_in = 1000;

    % not used
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;

    params('gamma') = 0.1;
    params('B') = 0.1;
    params('init_tau') = 1;
    params('init_alpha') = 1;

    [u, ~] = mcmc_gamma(params);
    u_avg = mean(sign(u(:, burn_in:end)), 2); %avg the rows

    %figure(1)
    %clf
    %scatter_twomoons_classify(data, u_avg, params('label_data'))
    p = count_correct(u_avg, params('label_data'), [zeros(floor(N/2)+1,1) - 1; zeros(N-(floor(N/2)+1),1) + 1]);
    %u_avg_accept = mean(u_accept(burn_in:end))

end