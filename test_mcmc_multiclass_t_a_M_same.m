function correct_p = test_mcmc_multiclass_t_a_M_same(percent_fidelity)

    params = containers.Map;
    params('laplacian') = "self tuning";
    
    digs = [1,3,4,5,9];
    saved = true;
    
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
    
    % not used
    %{
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;
    %}
    %
    
    params('num_iterations') = 100001;
    params('burn_in') = 2000;
    params('movie') = true;
    params('movie_often') = 2000; 
    
    params('gamma') = 0.0001;
    params('B') = 0.1;
    
    params('init_tau') = 0;
    params('tau_epsilon') = 0;
    params('tau_min') = 0;
    params('tau_max') = 20;
    
    
    params('init_alpha') = 1;
    params('alpha_epsilon') = 0;
    params('alpha_min') = 0.1;
    params('alpha_max') = 90;
    
    params('init_M') = 50;
    params('M_max_jump') = 0;
    params('M_min') = 1;
    params('M_max') = 50;
    
    params('remove-zero-eig') = true;
    
    correct_p = mcmc_multiclass_t_a_M_same(params);    
end
