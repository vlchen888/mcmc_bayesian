function cont = test_mcmc_hier_v_M(percent_fidelity, sigma)

    params = containers.Map;

    params('data_set') = "mnist";
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
        load('mnistdata.mat')
        params('data') = data;
        load('mnisttruth.mat')
        truth = truth(:, 1) - truth(:, 2);
        params('truth') = truth;
        params('label_data') = generate_fidelity(percent_fidelity, truth, length(data));
    end
    
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
    params('a') = 0.8;
    
    %params('epsilon') = 1e-13;
    
    params('epsilon') = .01;
    params('tau') = 1;
    params('alpha') = 35;
    
    params('init_M') = 50;
    params('max_M_jump') = 0;
    params('min_M') = 1;
    params('max_M') = 50;
        
    cont = mcmc_learn_v_M(params);
end