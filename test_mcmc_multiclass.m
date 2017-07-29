function p = test_mcmc_multiclass(percent_fidelity)

    params = containers.Map;

    params('laplacian') = "self tuning";
    
    digs = [4, 9];
    params('digs') = digs;
    k = 2;
    params('k') = k;
    load('mnist49data.mat')
    params('data') = data;
    load('mnist49truth.mat')
    params('truth') = truth';
    params('label_data') = generate_fidelity_multiclass(percent_fidelity, params('truth'), length(data), k);
    
    params('num_iterations') = 100000;
    burn_in = 5000;
    params('burn_in') = burn_in;

    % not used
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;
    %
    
    params('gamma') = 0.1;
    params('B') = 0.4;
    params('tau') = 2;
    params('alpha') = 35;
    
    
        
    [signs_all] ...
        = mcmc_multiclass(params);
    final_class = compute_S(mean((signs_all(:, :, burn_in:end)), 3));
    
    p = count_correct_multiclass(final_class, params('label_data'), params('truth'));
    
end

function S = compute_S(u)
    [~, I] = max(u,[],2);
    A = eye(length(u));
    S = A(:, I);
end