
function print_results_xi()
    
    params = containers.Map;
    params('data_set') = string('voting');
    
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
        load('moondata.mat')
        %%%% Currently need to set neg, pos manually.
        %%%% Should automate in the future.
        params('data') = data;
        set_pos     = 50:20:450;
        set_neg     = 550:20:950;
        params('set_neg') = set_neg;
        params('set_pos') = set_pos;
    end
    
    
    [num_data, ~] = size(params('data'));
    
    num_iterations = 1000;
    burn_in = 1;
    movie = 0;
    
    params('num_iterations') = num_iterations;
    params('burn_in') = burn_in;
    
    label_data = init(num_data, set_neg, set_pos);
    params('label_data') = label_data;
    
    params('p') = 2;
    params('q') = 2;
    params('l') = 1;

    params('gamma') = 0.0001;
    params('B')     = 0.1;
    
    params('init_tau')      = 20;
    params('init_alpha')    = 5;
    
    params('min_tau')      = 0.1;
    params('max_tau')    = 60;
    

    params('min_alpha')      = 0.1;
    params('max_alpha')    = 60;
    
    params('alpha_epsilon')      = 1;
    params('tau_epsilon')    = 1;
    
    start_time = cputime;
    [tau_all, alpha_all, std, xi_accept, tau_accept, alpha_accept] =...
        mcmc_learn_t_a_noncentered(params);
    elapsed_time = cputime - start_time;
    params('elapsed_time') = elapsed_time;
    
    %%%%% Take averages of tau, alpha over time as well %%%%%
    u_avg = zeros(num_data, num_iterations);
    tau_avg = zeros(1, num_iterations);
    alpha_avg = zeros(1, num_iterations);
    
    xi_accept_avg = zeros(1, num_iterations);
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
            
            xi_accept_avg(i+1) = ((i-burn_in)*xi_accept_avg(i) + xi_accept(i+1))/(i-burn_in+1);
            tau_accept_avg(i+1) = ((i-burn_in)*tau_accept_avg(i) + tau_accept(i+1))/(i-burn_in+1);
            alpha_accept_avg(i+1) = ((i-burn_in)*alpha_accept_avg(i) + alpha_accept(i+1))/(i-burn_in+1);
        end
        
        %%%%% CODE FOR AVG MOVIE %%%%
        if movie
            often = floor(num_iterations/100);
            if mod(i, often) == 1
                %%% Plot trace of u? %%%
                clf
                subplotBar(std(:, i))

                subplot(2,2,2)
                subplot_scatter_twomoons_classify(data, u_avg(:,i), label_data);

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
    params('xi_accept_avg') = xi_accept_avg;
    params('tau_accept_avg') = tau_accept_avg;
    params('alpha_accept_avg') = alpha_accept_avg;
    print_figures(params)
    
    params('tau_avg') = tau_avg(num_iterations);
    params('alpha_avg') = alpha_avg(num_iterations);
    
    if params('data_set') == string('voting')
        params('correct_percent') = count_correct(u_avg(:, num_iterations), set_neg, set_pos, [zeros(267,1) - 1; zeros(168,1) + 1]);
    elseif params('data_set') == string('moon')
        params('correct_percent') = count_correct(u_avg(:, num_iterations), set_neg, set_pos, ...
        [zeros(num_data/2,1) + 1; zeros(num_data/2, 1) - 1]);
    end
    
    params('run_description') = 'MCMC self-tuning Laplacian, xi tau alpha parameterization.\n';
    
    print_info_file(params);
end

function print_figures(params)
    num_iterations = params('num_iterations');
    tau_all = params('tau_all');
    alpha_all = params('alpha_all');
    u_avg = params('u_avg');
    tau_avg = params('tau_avg');
    alpha_avg = params('alpha_avg');
    xi_accept_avg = params('xi_accept_avg');
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
    plotBar(u_avg(:, num_iterations), 1)
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
    plot(burn_in+1:num_iterations, xi_accept_avg(burn_in+1:end))
    ylabel('\xi acceptance probability');
    fname = 'print_runs/acceptance_xi_probability.png';
    print('-r144','-dpng',fname);
    
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

function subplot_scatter_twomoons_classify(data, final_avg, label_data)
    colormap(redbluecmap(5))
    colors = -sign(final_avg); % to make blue = +, red = -
    scatter(data(:,1), data(:,2), 5 , colors)
    hold on
    for i = 1:length(label_data)
        if label_data(i) == 1
            scatter(data(i,1),data(i,2), 20, 'b', 'd', 'filled');
        elseif label_data(i) == -1
            scatter(data(i,1),data(i,2), 20, 'r', 'd', 'filled');
        end
    end
    hold off
end

function scatter_twomoons_classify(data, final_avg, label_data)
    colormap(redbluecmap(5))
    colors = -sign(final_avg); % to make blue = +, red = -
    scatter(data(:,1), data(:,2), 5 , colors)
    hold on
    for i = 1:length(label_data)
        if label_data(i) == 1
            scatter(data(i,1),data(i,2), 20, 'b', 'd', 'filled');
        elseif label_data(i) == -1
            scatter(data(i,1),data(i,2), 20, 'r', 'd', 'filled');
        end
    end
    hold off
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
elapsed_time = params('elapsed_time');
init_tau = params('init_tau');
init_alpha = params('init_alpha');

fileID = fopen('print_runs/run_info.txt','w');
fprintf(fileID, 'This is an autogenerated text file with run info and summary!\n');
fprintf(fileID, ['Description: ' run_description]);
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
fprintf(fileID, 'Time elapsed: %.2f s\n', elapsed_time);

end

function subplotBar(avg)
    num_pos = 0;
    for i=1:length(avg)
        if avg(i) > 0
            num_pos = num_pos + 1;
        end
    end

    posx = zeros(num_pos,1);
    posy = zeros(num_pos,1);
    ppos = 1;
    negx = zeros(435-num_pos,1);
    negy = zeros(435-num_pos,1);
    pneg = 1;
    for i=1:length(avg)
        if avg(i) > 0
            posx(ppos) = i;
            posy(ppos) = avg(i);
            ppos = ppos + 1;
        else
            negx(pneg) = i;
            negy(pneg) = avg(i);
            pneg = pneg + 1;
        end
    end
    
    
    subplot(2,2,1)
    hold on
    bar(posx, posy, 'b')
    bar(negx, negy, 'r')
    hold off
end

function plotBar(avg, k)
    
    num_pos = 0;
    for i=1:length(avg)
        if avg(i) > 0
            num_pos = num_pos + 1;
        end
    end

    posx = zeros(num_pos,1);
    posy = zeros(num_pos,1);
    ppos = 1;
    negx = zeros(435-num_pos,1);
    negy = zeros(435-num_pos,1);
    pneg = 1;
    for i=1:length(avg)
        if avg(i) > 0
            posx(ppos) = i;
            posy(ppos) = avg(i);
            ppos = ppos + 1;
        else
            negx(pneg) = i;
            negy(pneg) = avg(i);
            pneg = pneg + 1;
        end
    end
    
    
    figure(k)
    clf
    hold on
    bar(posx, posy, 'b')
    bar(negx, negy, 'r')
    hold off
end
