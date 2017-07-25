
function [tau_all, alpha_all] = print_results()
%For printing traces, acceptance probabilities, run statistics, etc.
%   Set run parameters here.
    params = containers.Map;
    params('data_set') = string('voting');
    
    params('algorithm') = string('noncentered');
    params('laplacian') = string('un');
    
    num_iterations = 100000;
    burn_in = 10000;
    params('num_iterations') = num_iterations;
    params('burn_in') = burn_in;
    
    movie = 0;
    
    params('p') = 2;
    params('q') = 2;
    params('l') = 1;

    %%%% CENTERED PARAMS for voting %%%%
    %{
    params('gamma')         = 0.0001;
    params('B')             = 0.4;
    
    params('init_tau')      = 20;
    params('init_alpha')    = 20;
    
    params('min_tau')       = 0.1;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;
    
    params('alpha_epsilon') = 0.1;
    params('tau_epsilon')   = 1;
    %}
    
    %%%% CENTERED PARAMS for intertwined moons %%%%
    %{
    params('gamma')         = 0.1;
    params('B')             = 0.4;
    params('init_tau')      = 20;
    params('init_alpha')    = 20;

    params('min_tau')       = 0.01;
    params('max_tau')       = 60;

    params('min_alpha')     = 1;
    params('max_alpha')     = 60;

    params('alpha_epsilon') = 0.5;
    params('tau_epsilon')   = 1;
    %}
        
    %%%% NONCENTEREED PARAMS for voting %%%%
    
    params('gamma')         = 0.0001;
    params('B')             = 0.1;

    params('init_tau')      = 20;
    params('init_alpha')    = 20;

    params('min_tau')       = 0.1;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;

    params('alpha_epsilon') = 1;
    params('tau_epsilon')   = 1;
    
    
    %%%% NONCENTERED PARAMS for moons %%%%
    %{
    params('gamma')         = 0.1;
    params('B')             = 0.4;
    params('init_tau')      = 20;
    params('init_alpha')    = 20;

    params('min_tau')       = 0.01;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;

    params('alpha_epsilon') = 0.5;
    params('tau_epsilon')   = 1;
    %}
    
    %%%% NONCENTERED PARAMS intertwined moons, self tuning %%%%
    %%%% Gets convergence in ~400 iterations, ~98% accuracy on sigma = 0.1
    %{
    params('gamma')         = 0.1;
    params('B')             = 0.4;
    
    params('init_tau')      = 1;
    params('init_alpha')    = 1;
    
    params('min_tau')       = 0.01;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;
    
    params('alpha_epsilon') = .5;
    params('tau_epsilon')   = .3;
    %}
    
    %%%% NONCENTERED PARAMS for intertwined moons, unnormalized%%%%
    %{
    params('gamma')         = 0.0001;
    params('B')             = 0.002;
    
    params('init_tau')      = 20;
    params('init_alpha')    = 20;
    
    params('min_tau')       = 0.1;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;
    
    params('alpha_epsilon') = 0.1;
    params('tau_epsilon')   = 1;
    %}
    
    %%%% t,a,M PARAMS for intertwined moons %%%%
    %{
    params('gamma')         = 0.0001;
    params('B')             = 0.1;
    
    params('init_tau')      = 20;
    params('init_alpha')    = 20;
    params('init_M')        = 500;
    
    params('min_tau')       = 0.1;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;
    
    params('min_M')       = 1;
    params('max_M')       = 1000;
    
    params('alpha_epsilon') = 0.1;
    params('tau_epsilon')   = 0.1;
    %}
    
    if params('data_set') == string('voting')
        load('data3.mat')
        data = X;
        params('data') = data;
        params('truth') = [-ones(267,1); ones(168,1)];
        %percent_fidelity = .0115;
        percent_fidelity = .08;
        params('label_data') = generate_fidelity(percent_fidelity, params('truth'), length(data));        
    elseif params('data_set') == string('moons')
        load('intertwine_moon.mat')
        %%%% Currently need to set neg, pos manually.
        %%%% Should automate in the future.
        data = d;
        params('data') = data;
        params('truth') = [-ones(floor(N/2)+1,1); ones(N-(floor(N/2)+1),1)];
        params('label_data') = generate_fidelity(percent_fidelity, params('truth'), length(data));        
    end    
    
    if params('algorithm') == string('noncentered')
        [tau_all, alpha_all, std, var_accept, tau_accept, alpha_accept] =...
            mcmc_learn_t_a_noncentered(params);
    elseif params('algorithm') == string('centered')
        [tau_all, alpha_all, std, var_accept, tau_accept, alpha_accept] =...
            mcmc_learn_t_a(params);
    end
    
    %%%%% Take averages of tau, alpha over time as well %%%%%    
    u_avg = movmean(sign(std(:,burn_in:end)), [num_iterations-burn_in+1,0], 2);
    u_var = 1 - u_avg.^2;
    tau_avg = movmean(tau_all(burn_in:end), [num_iterations-burn_in+1,0]);
    alpha_avg = movmean(alpha_all(burn_in:end), [num_iterations-burn_in+1,0]);
    var_accept_avg = movmean(var_accept, [num_iterations,0]);
    tau_accept_avg = movmean(tau_accept, [num_iterations,0]);
    alpha_accept_avg = movmean(alpha_accept, [num_iterations,0]);
    
    if movie
        for i=1:floor(num_iterations/100):num_iterations-1
            %%%%% CODE FOR AVG MOVIE %%%%
            clf
            %%% Plot trace of u? %%%
            subplot(2,2,1)
            plotBar(std(:, i))

            %subplotBar(u_avg(:, i))

            if params('data_set') == string('moons')
                subplot(2,2,2)
                scatter_twomoons_classify(data, u_avg(:,i), params('label_data'));
            end

            subplot(2,2,3)
            plot(1:i, tau_all(1:i));
            ylabel('\tau trace');

            subplot(2,2,4)
            plot(1:i, alpha_all(1:i));
            ylabel('\alpha trace');

            fname = sprintf('figs/step_%i.png',floor(i/often) + 1);
            print('-r144','-dpng',fname);
        end
    end
    params('tau_all') = tau_all;
    params('alpha_all') = alpha_all;
    params('u_avg') = u_avg;
    params('tau_avg') = tau_avg;
    params('alpha_avg') = alpha_avg;
    params('var_accept_avg') = var_accept_avg;
    params('tau_accept_avg') = tau_accept_avg;
    params('alpha_accept_avg') = alpha_accept_avg;
    print_figures(params)
    
    params('tau_avg') = tau_avg(end);
    params('alpha_avg') = alpha_avg(end);
    
    params('correct_percent') = count_correct(u_avg(:, end), params('label_data'), params('truth'));
    
    if params('algorithm') == string('noncentered')
        params('run_description') = 'MCMC self-tuning Laplacian, noncentered parameterization.\n';
    elseif params('algorithm') == string('centered')
        params('run_description') = 'MCMC self-tuning Laplacian, centered parameterization.\n';
    end
        
    print_info_file(params);
end

function print_figures(params)
    num_iterations = params('num_iterations');
    tau_all = params('tau_all');
    alpha_all = params('alpha_all');
    u_avg = params('u_avg');
    tau_avg = params('tau_avg');
    alpha_avg = params('alpha_avg');
    var_accept_avg = params('var_accept_avg');
    tau_accept_avg = params('tau_accept_avg');
    alpha_accept_avg = params('alpha_accept_avg');
    burn_in = params('burn_in');
    data = params('data');
    
    figure(1)
    set(gcf, 'Position', [100, 300, 800, 300])
    
    %%%%% PRINT TRACES OF TAU AND ALPHA %%%%%
    clf
    plot(1:num_iterations, tau_all, burn_in:num_iterations, tau_avg);
    legend('\tau trace', '\tau running average');
    fname = 'print_runs/trace_tau.png';
    print('-r144','-dpng',fname);
    
    clf
    plot(1:num_iterations, alpha_all,burn_in:num_iterations, alpha_avg);
    legend('\alpha trace', '\alpha running average');
    fname = 'print_runs/trace_alpha.png';
    print('-r144','-dpng',fname);
    
    %%%%%%% FINAL AVG %%%%%%%
    clf
    plotBar(u_avg(:, end))
    fname = 'print_runs/final_avg.png';
    print('-r144','-dpng',fname);
    
    %%%%%%% TAU AVG %%%%%%%
    %{
    clf
    plot(burn_in+1:num_iterations, tau_avg(burn_in+1:num_iterations))
    ylabel('\tau running average');
    fname = 'print_runs/avg_tau.png';
    print('-r144','-dpng',fname);
    
    %%%%%%% ALPHA AVG %%%%%%%
    clf
    plot(burn_in+1:num_iterations, alpha_avg(burn_in+1:num_iterations))
    ylabel('\alpha running average');
    fname = 'print_runs/avg_alpha.png';
    print('-r144','-dpng',fname);
    %}
    
    clf
    if params('algorithm') == string('noncentered')
        plot(var_accept_avg)
        ylabel('\xi acceptance probability');
        fname = 'print_runs/acceptance_xi_probability.png';
        print('-r144','-dpng',fname);
    elseif params('algorithm') == string('centered')
        plot(var_accept_avg)
        ylabel('u acceptance probability');
        fname = 'print_runs/acceptance_u_probability.png';
        print('-r144','-dpng',fname);
    end
    
    clf
    plot(tau_accept_avg)
    ylabel('\tau acceptance probability');
    fname = 'print_runs/acceptance_tau_probability.png';
    print('-r144','-dpng',fname);
    
    clf
    plot(alpha_accept_avg)
    ylabel('\alpha acceptance probability');
    fname = 'print_runs/acceptance_alpha_probability.png';
    print('-r144','-dpng',fname);
    
    if params('data_set') == string('moons')
        clf
        scatter_twomoons_classify(data, u_avg(:, end), params('label_data'))
        fname = 'print_runs/final_scatter.png';
        print('-r144','-dpng',fname);
    end
end

function print_info_file(params)
run_description = params('run_description');
num_iterations = params('num_iterations');
burn_in = params('burn_in');
p = params('p');
q = params('q');
l = params('l');
B = params('B');
gamma = params('gamma');
tau_avg = params('tau_avg');
alpha_avg = params('alpha_avg');
label_data = params('label_data');
correct_percent = params('correct_percent');
tau_epsilon = params('tau_epsilon');
alpha_epsilon = params('alpha_epsilon');
init_tau = params('init_tau');
init_alpha = params('init_alpha');

fileID = fopen('print_runs/run_info.txt','w');
fprintf(fileID, 'This is an autogenerated text file with run info and summary!\n');
fprintf(fileID, ['Description:\n ' run_description]);
fprintf(fileID, 'Laplacian: %s\n', params('laplacian'));
fprintf(fileID, 'Iterations = %d, Burn in period = %d\n', num_iterations, burn_in);
fprintf(fileID, 'Parameters of weight function: p = %d, q = %d, l = %d\n', p, q, l);
fprintf(fileID, 'Beta = %f, Gamma = %d\n', B, gamma);
fprintf(fileID, 'Tau epsilon = %f, Alpha epsilon = %f\n', tau_epsilon, alpha_epsilon);
fprintf(fileID, 'Initial tau = %f, Initial alpha = %f\n', init_tau, init_alpha);
fprintf(fileID, 'Average tau = %f, Average alpha = %f\n', tau_avg, alpha_avg);
fprintf(fileID, 'Label Data:\n');
fprintf(fileID, '+1: %d\n', find(label_data>0));
fprintf(fileID, '-1: %d\n', find(label_data<0));
fprintf(fileID, 'Percent correctly classified: %f\n', correct_percent);

end