function p = test_mcmc_gamma_accuracy()
rng(1)
params = containers.Map;
params('data_set') = string('moons');
params('laplacian') = string('self tuning');
params('percent_fidelity') = 0.005;

    
params('num_iterations') = 100000;
burn_in = 1000;

if params('data_set') == string('voting')
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;
    params('gamma') = 0.0001;
    params('B') = 0.6;
    params('init_tau') = 2.4;
    params('init_alpha') = 60;
    
    load('data3.mat')
    data = X;
    params('data') = data;
    params('label_data') = generate_voting_fidelity(params('percent_fidelity'));

elseif params('data_set') == string('moons')
    % not used
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;
    
    params('gamma') = 0.1;
    params('B') = 0.4;
    params('init_tau') = 0.1;
    params('init_alpha') = 1;

    N = 2000;
    sigma = 0.02;
    data = moondata(1,100,N,sigma);
    params('data') = data;
    params('label_data') = generate_moons_fidelity(params('percent_fidelity'), N);
end

start_time = cputime;
[u, u_accept] = mcmc_gamma(params);
elapsed_time = cputime - start_time
u_avg = mean(sign(u(:, burn_in:end)), 2); %avg the rows

if params('data_set') == string('voting')
    figure(1)
    clf
    plotBar(u_avg);
    p = count_correct(u_avg, params('label_data'), [zeros(267,1) - 1; zeros(168,1) + 1]);
elseif params('data_set') == string('moons')
    figure(1)
    clf
    scatter_twomoons_classify(data, u_avg, params('label_data'))
    p = count_correct(u_avg, params('label_data'), [zeros(floor(N/2)+1,1) - 1; zeros(N-(floor(N/2)+1),1) + 1]);
end
u_avg_accept = mean(u_accept(burn_in:end))
end