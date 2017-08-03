function p = test_mcmc_multiclass_t_a(percent_fidelity)

    params = containers.Map;

    params('laplacian') = "self tuning";
    
    digs = [1, 4, 7, 9];
    saved = false;
    
    params('digs') = digs;
    k = length(digs);
    params('k') = k;
    if saved
        load('mnistdata.mat')
        params('data') = data;
        load('mnisttruth.mat')
        params('truth') = truth;
    else
        [data, truth] = generate_mnist_data(digs);
        save('mnistdata.mat','data');
        save('mnisttruth.mat','truth');
        params('data') = data;
        params('truth') = truth;
    end
    params('label_data') = generate_fidelity_multiclass(percent_fidelity, params('truth'), length(data), k);
    
    params('num_iterations') = 100001;
    burn_in = 5000;
    params('burn_in') = burn_in;

    % not used
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;
    %
    
    params('gamma') = 0.1;
    params('B') = 0.1;
    params('tau') = 2;
    params('alpha') = 35;
    
    params('init_tau') = 2;
    params('init_alpha') = 20;
    
    params('tau_epsilon') = 0.1;
    params('alpha_epsilon') = 0.5;
    
    params('tau_min') = 0.1;
    params('tau_max') = 50;
    params('alpha_min') = 0.1;
    params('alpha_max') = 60;
    
    [u_all] = mcmc_multiclass_t_a(params);
    final_class = compute_S_multiclass(mean((u_all(:, :, burn_in:end)), 3), k);
    
    p = count_correct_multiclass(final_class, params('label_data'), params('truth'));
    
end