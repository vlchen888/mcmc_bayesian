function p = test_mcmc_multiclass(percent_fidelity)

    params = containers.Map;

    params('laplacian') = "self tuning";
    
    digs = [1, 4, 9];
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
    
    u_avg = mcmc_multiclass(params);
    
    p = count_correct_multiclass(compute_S_multiclass(u_avg,k), params('label_data'), params('truth'));
    
end