function cont = test_multiclass_t_a_M(percent_fidelity)
% For parameter settings and testing multiclass MCMC

    params = containers.Map;
    params('laplacian') = "self tuning";
    
    digs = [3,4,5,9];
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
    
    params('num_iterations') = 25001;
    params('burn_in') = 2000;
    params('movie') = false;
    params('movie_period') = 25000;
    params('accuracy_period') = 5000;
    
    params('gamma') = .0001;
   
    params('B_init') = 0.1;
    params('B_update_period') = 500;
    params('B_target_p') = 0.5;
    params('B_burn_in') = 10000;
    
    
    params('init_tau') = 2;
    params('tau_epsilon') = 0.2;
    params('tau_min') = -6;
    params('tau_max') = 6;
    
    
    params('init_alpha') = 1;
    params('alpha_epsilon') = 5;
    params('alpha_min') = 0.1;
    params('alpha_max') = 110;
    
    params('init_M') = 50;
    params('M_max_jump') = 15;
    params('M_min') = 1;
    params('M_max') = 80;
    
    params('remove-zero-eig') = true;
    params('adaptive') = true;
    
    cont = mcmc_multiclass_t_a_M_same(params);    
end
