function correct_p = test_mcmc_multiclass_t_a_M_same(percent_fidelity)

    params = containers.Map;
    params('laplacian') = "self tuning";
    
    digs = [1,3,4,5,9];
    saved = true;
    
    params('digs') = digs;
    params('movie') = false;
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
    
    params('num_iterations') = 25001;
    burn_in = 5000;
    params('burn_in') = burn_in;

    % not used
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;
    %
    
    params('gamma') = 0.0001;
    params('B') = 0.1;
    
    params('init_tau') = 0.7;
    params('tau_epsilon') = 0;
    params('tau_min') = 0.01;
    params('tau_max') = 20;
    
    
    params('init_alpha') = 40;
    params('alpha_epsilon') = 1.3;
    params('alpha_min') = 0.1;
    params('alpha_max') = 90;
    
    params('init_M') = 50;
    params('M_max_jump') = 20;
    params('M_min') = 1;
    params('M_max') = 70;
    
    correct_p = mcmc_multiclass_t_a_M_same(params);    
end
