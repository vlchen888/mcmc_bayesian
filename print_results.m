
function [tau_all, alpha_all] = print_results()
%For printing traces, acceptance probabilities, run statistics, etc.
%   Set run parameters here.
    params = containers.Map;
    params('data_set') = string('moons');
    
    params('algorithm') = string('t a M');
    params('laplacian') = string('self tuning');
    
    num_iterations = 10000;
    burn_in = 1000;
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
    params('gamma')         = 0.0001;
    params('B')             = 0.1;
    
    params('init_tau')      = 20;
    params('init_alpha')    = 20;
    
    params('min_tau')       = 0.1;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;
    
    params('alpha_epsilon') = 0.1;
    params('tau_epsilon')   = 0.1;
    %}
        
    %%%% NONCENTEREED PARAMS for voting %%%%
    %{
    params('gamma')         = 0.0001;
    params('B')             = 0.1;

    params('init_tau')      = 30;
    params('init_alpha')    = 20;

    params('min_tau')       = 0.1;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;

    params('alpha_epsilon') = 1;
    params('tau_epsilon')   = 1;
    %}
    
    %%%% NONCENTERED PARAMS for moons %%%%
    %{
    params('gamma')         = 2;
    params('B')             = 0.0003;
    
    params('init_tau')      = 40;
    params('init_alpha')    = 20;
    
    params('min_tau')       = 0.1;
    params('max_tau')       = 60;

    params('min_alpha')     = 0.1;
    params('max_alpha')     = 60;
    
    params('alpha_epsilon') = 2;
    params('tau_epsilon')   = 3;
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
    
    
    if params('data_set') == string('voting')
        load('data3.mat')
        data = X;
        params('data') = data;
        % set_neg     = [25 26 27 28];
        % set_pos     = [280 281 282 283];
        set_neg     = 20:30;
        set_pos     = 280:290;
        params('set_neg') = set_neg;
        params('set_pos') = set_pos;
        
    elseif params('data_set') == string('moons')
        load('intertwine_moon.mat')
        %%%% Currently need to set neg, pos manually.
        %%%% Should automate in the future.
        data = d;
        params('data') = data;
        set_pos     = 50:20:450;
        set_neg     = 550:20:950;
        params('set_neg') = set_neg;
        params('set_pos') = set_pos;
    end
    
    [num_data, ~] = size(params('data'));
    
    label_data = init(num_data, set_neg, set_pos);
    params('label_data') = label_data;
    
    if params('algorithm') == string('noncentered')
        [tau_all, alpha_all, std, var_accept, tau_accept, alpha_accept] =...
            mcmc_learn_t_a_noncentered(params);
    elseif params('algorithm') == string('centered')
        [tau_all, alpha_all, std, var_accept, tau_accept, alpha_accept] =...
            mcmc_learn_t_a(params);
    end
    
    %%%%% Take averages of tau, alpha over time as well %%%%%
    u_avg = zeros(num_data, num_iterations);
    tau_avg = zeros(1, num_iterations);
    alpha_avg = zeros(1, num_iterations);
    
    var_accept_avg = zeros(1, num_iterations);
    tau_accept_avg = zeros(1, num_iterations);
    alpha_accept_avg = zeros(1, num_iterations);
    
    %%%%%%% Do running average and variance as well %%%%%%
    % u_avg(j, n) represents the average of the sign of the first n clusterings
    u_var = zeros(num_data, num_iterations);
 
    for i=1:num_iterations-1
        if i >= burn_in
            u_avg(:, i+1) = ((i-burn_in)*u_avg(:, i) + sign(std(:, i+1)))/(i-burn_in+1);
            u_var(:, i+1) = 1-u_avg(:, i+1).^2;
            
            tau_avg(i+1) = ((i-burn_in)*tau_avg(i) + tau_all(i+1))/(i-burn_in+1);
            alpha_avg(i+1) = ((i-burn_in)*alpha_avg(i) + alpha_all(i+1))/(i-burn_in+1);
            
            var_accept_avg(i+1) = ((i-burn_in)*var_accept_avg(i) + var_accept(i+1))/(i-burn_in+1);
            tau_accept_avg(i+1) = ((i-burn_in)*tau_accept_avg(i) + tau_accept(i+1))/(i-burn_in+1);
            alpha_accept_avg(i+1) = ((i-burn_in)*alpha_accept_avg(i) + alpha_accept(i+1))/(i-burn_in+1);
        end
        
        %%%%% CODE FOR AVG MOVIE %%%%
        if movie
            often = floor(num_iterations/100);
            if mod(i-1, often) == 0
                
                clf
                %%% Plot trace of u? %%%
                subplot(2,2,1)
                plotBar(std(:, i))
                
                %subplotBar(u_avg(:, i))
                
                if params('data_set') == string('moons')
                    subplot(2,2,2)
                    scatter_twomoons_classify(data, u_avg(:,i), label_data);
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
    
    params('tau_avg') = tau_avg(num_iterations);
    params('alpha_avg') = alpha_avg(num_iterations);
    
    if params('data_set') == string('voting')
        params('correct_percent') = count_correct(u_avg(:, num_iterations), label_data, [zeros(267,1) - 1; zeros(168,1) + 1]);
    elseif params('data_set') == string('moons')
        params('correct_percent') = count_correct(u_avg(:, num_iterations), label_data, ...
        [zeros(num_data/2,1) + 1; zeros(num_data/2, 1) - 1]);
    end
    
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
    label_data = params('label_data');
    
    %%%%% PRINT TRACES OF TAU AND ALPHA %%%%%
    clf
    plot(1:num_iterations, tau_all);
    ylabel('\tau trace');
    fname = 'print_runs/trace_tau.png';
    print('-r144','-dpng',fname);
    
    clf
    plot(1:num_iterations, alpha_all);
    ylabel('\alpha trace');
    fname = 'print_runs/trace_alpha.png';
    print('-r144','-dpng',fname);
    
    %%%%%%% FINAL AVG %%%%%%%
    clf
    plotBar(u_avg(:, num_iterations))
    fname = 'print_runs/final_avg.png';
    print('-r144','-dpng',fname);
    
    %%%%%%% TAU AVG %%%%%%%
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
    
    clf
    if params('algorithm') == string('noncentered')
        plot(burn_in+1:num_iterations, var_accept_avg(burn_in+1:end))
        ylabel('\xi acceptance probability');
        fname = 'print_runs/acceptance_xi_probability.png';
        print('-r144','-dpng',fname);
    elseif params('algorithm') == string('centered')
        plot(burn_in+1:num_iterations, var_accept_avg(burn_in+1:end))
        ylabel('u acceptance probability');
        fname = 'print_runs/acceptance_u_probability.png';
        print('-r144','-dpng',fname);
    end
    
    clf
    plot(burn_in+1:num_iterations, tau_accept_avg(burn_in+1:end))
    ylabel('\tau acceptance probability');
    fname = 'print_runs/acceptance_tau_probability.png';
    print('-r144','-dpng',fname);
    
    clf
    plot(burn_in+1:num_iterations, alpha_accept_avg(burn_in+1:end))
    ylabel('\alpha acceptance probability');
    fname = 'print_runs/acceptance_alpha_probability.png';
    print('-r144','-dpng',fname);
    
    if params('data_set') == string('moons')
        clf
        scatter_twomoons_classify(data, u_avg(:, num_iterations), label_data)
        fname = 'print_runs/final_scatter.png';
        print('-r144','-dpng',fname);
    end
end

function u = init(num_senators, set_neg, set_pos)
    u = zeros(num_senators, 1);
    
    % label some of the data
    for i=1:length(set_pos)
        u(set_pos(i))=1;
    end
    for i=1:length(set_neg)
        u(set_neg(i))=-1;
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
set_neg = params('set_neg');
set_pos = params('set_pos');
correct_percent = params('correct_percent');
tau_epsilon = params('tau_epsilon');
alpha_epsilon = params('alpha_epsilon');
init_tau = params('init_tau');
init_alpha = params('init_alpha');

fileID = fopen('print_runs/run_info.txt','w');
fprintf(fileID, 'This is an autogenerated text file with run info and summary!\n');
fprintf(fileID, ['Description: ' run_description]);
fprintf(fileID, 'Laplacian: %s', params('laplacian'));
fprintf(fileID, 'Iterations = %d, Burn in period = %d\n', num_iterations, burn_in);
fprintf(fileID, 'Parameters of weight function: p = %d, q = %d, l = %d\n', p, q, l);
fprintf(fileID, 'Beta = %f, Gamma = %d\n', B, gamma);
fprintf(fileID, 'Tau epsilon = %f, Alpha epsilon = %f\n', tau_epsilon, alpha_epsilon);
fprintf(fileID, 'Initial tau = %f, Initial alpha = %f\n', init_tau, init_alpha);
fprintf(fileID, 'Average tau = %f, Average alpha = %f\n', tau_avg, alpha_avg);
fprintf(fileID, 'Label Data:\n');
fprintf(fileID, '+1: %d\n', set_pos);
fprintf(fileID, '-1: %d\n', set_neg);
fprintf(fileID, 'Percent correctly classified: %f\n', correct_percent);

end