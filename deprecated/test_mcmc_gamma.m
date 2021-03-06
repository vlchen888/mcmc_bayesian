function p = test_mcmc_gamma(percent_fidelity, sigma)

    params = containers.Map;
    params('data_set') = "moons";
    params('laplacian') = "self tuning";

    if params('data_set') == "moons"
        
        N = 2000;
        data = moondata(1,100,N,sigma);
        params('truth') = [-ones(floor(N/2)+1,1); ones(N-(floor(N/2)+1),1)];
        params('data') = data;
        params('label_data') = generate_fidelity(percent_fidelity, params('truth'), length(data));
        
        %{
        N = 1000;
        load('intertwine_moon.mat')
        data = d;
        params('data') = data;
        set_pos     = 50:20:450;
        set_neg     = 550:20:950;
        params('label_data') = init(length(data),set_neg,set_pos);
        %}

    elseif params('data_set') == "voting"
        load('data3.mat')
        data = X;
        params('data') = data;
        params('truth') = [-ones(267,1); ones(168,1)];
        params('label_data') = generate_fidelity(percent_fidelity, params('truth'), length(data));
    end

    figure(3)
    plot(params('label_data'))
    
    params('num_iterations') = 100000;
    burn_in = 5000;

    % not used
    params('p') = 2;
    params('q') = 2;
    params('l') = 1.25;
    %

    
    
    params('gamma') = 0.1;
    params('B') = 0.4;
    params('init_tau') = 2;
    params('init_alpha') = 35;
    
    %{
    params('gamma')         = 0.0001;
    params('B')             = 0.4;
    
    params('init_tau')      = 1;
    params('init_alpha')    = 1;
    %}

    [u_all, u_accept] = mcmc_gamma(params);
    u_avg = mean(u_all(:, burn_in:end), 2); %avg the rows

    
    if params('data_set') == "moons"
        figure(1)
        set(gcf, 'Position', [100, 300, 800, 300])
        plotBar(u_avg)
        figure(2)
        set(gcf, 'Position', [100, 300, 800, 300])
        plot(movmean(u_accept,[length(u_accept) 0]))
        ylabel('u acceptance probability')
        figure(4)
        set(gcf, 'Position', [100, 300, 600, 500])
        scatter_twomoons_classify(data, u_avg, params('label_data'))
        
        p = count_correct(u_avg, params('label_data'), params('truth'));
    elseif params('data_set') == "voting"
        figure(1)
        set(gcf, 'Position', [100, 300, 800, 300])
        plotBar(u_avg)
        figure(2)
        set(gcf, 'Position', [100, 300, 800, 300])
        plot(movmean(u_accept,[length(u_accept) 0]))
        ylabel('u acceptance probability')
        figure(3)
        set(gcf, 'Position', [100, 300, 800, 300])
        plot_me = [1:5,271:275];
        plot_senator_traces(u_all, plot_me)
            
        p = count_correct(u_avg, params('label_data'), [-ones(267,1); ones(168,1)]);

        
    end
end