function p = test_mcmc_hier_v(percent_fidelity, sigma)

    params = containers.Map;

    params('data_set') = "moons";
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
        load('mnist49data.mat')
        params('data') = data;
        load('mnist49truth.mat')
        truth = truth(:, 1) - truth(:, 2);
        params('truth') = truth;
        params('label_data') = generate_fidelity(percent_fidelity, truth, length(data));
    end
    
    figure(1)
    plot(params('label_data'))

    params('num_iterations') = 100000;
    burn_in = 5000;
    params('burn_in') = burn_in;

    % not used
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;
    %
    
    params('gamma') = 0.1;
    params('B') = 0.1;
    params('a') = 0.5;
    params('epsilon') = 1e-13;
    params('tau') = 2;
    params('alpha') = 35;
    
    
        
    [std, v_all, xi_all, v_accept, xi_accept] ...
        = mcmc_learn_v(params);
    u_avg = mean(std(:, burn_in:end), 2); %avg the rows
    
    p = count_correct(u_avg, params('label_data'), params('truth'));
        
    figure(2)
    if params('data_set') == "moons"
        set(gcf, 'Position', [100, 300, 600, 500])
        scatter_twomoons_classify(data, u_avg, params('label_data'))
    elseif params('data_set') == "voting"
        set(gcf, 'Position', [100, 300, 800, 300])
        plotBar(u_avg)
    end
    
    figure(3)
    set(gcf, 'Position', [100, 300, 800, 300])
    plot(mean(xi_all.*v_all, 2));
    title('u_j average')
        
    figure(4)
    set(gcf, 'Position', [100, 300, 800, 300])
    subplot(211)
    plot(movmean(xi_accept, [length(xi_accept) 0]));
    title('\xi acceptance probability')
    subplot(212)
    plot(movmean(v_accept, [length(v_accept) 0]));
    title('v acceptance probability')

end