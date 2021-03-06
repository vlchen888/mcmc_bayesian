
function print_results_xi()
    
    %%%% Test voter data %%%%
    load('data3.mat')
    data = X;
    
    %%%% Test two moons data %%%%
    %d=100; N=1000; sigma=0.1;
    %data = moondata(d,N,sigma);
    %load('moondata.mat');
    %N = 1000;
    
    [num_data, ~] = size(data);
    num_iterations = 10000;
    burn_in = 500;
    
    %%%%  Voter data labeling %%%%
    set_neg     = [25 26 27 28];
    set_pos     = [280 281 282 283];
    
    % set_neg     = 20:30;
    % set_pos     = 280:290;
    
    %%%% Two moons labeling %%%% Currently need to set neg, pos manually.
    %%%% Should automate in the future.
    %set_pos     = 50:20:450;
    %set_neg     = 550:20:950;
    
    label_data = init(num_data, set_neg, set_pos);
    
    p = 2;
    q = 2;
    l = 1;

    gamma       = 2;
    B           = 0.0003;
    
    init_tau        = 30;
    init_alpha      = 5;
    
    min_tau         = 0.1;
    max_tau         = 60;
    
    min_alpha       = 0.1;
    max_alpha       = 60;
    
    alpha_epsilon   = 2;    % Jump alpha
    tau_epsilon     = 3;    % Jump tau
    
    start_time = cputime;
    [tau_all, alpha_all, std, xi_accept, tau_accept, alpha_accept] =...
        mcmc_learn_t_a_noncentered(data, num_iterations, label_data, p, q, l, ...
        gamma, B, init_tau, init_alpha, min_tau, max_tau, min_alpha, ...
        max_alpha, alpha_epsilon, tau_epsilon);
    elapsed_time = cputime - start_time;

    
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
        often = floor(num_iterations/100);
        if mod(i, often) == 1
            %%% Plot trace of u? %%%
            clf
            subplotBar(std(:, i))
           
            %subplot(2,2,2)
            %subplot_scatter_twomoons_classify(data, u_avg(:,i), label_data);
           
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
    
    %plot_me = [(1:10)'; (268:277)'];

    %%%%% PRINT RUNNING AVERAGE OF U %%%%%
    %plotAvg(f, num_iterations, plot_me, 1)
    %fname = 'print_runs/running_avg.png';
    %print('-r144','-dpng',fname);
    
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
    
    %%%%%%% Plot Tau? %%%%%%
    %figure(3)
    %histogram(tau_all)
    
    %figure(4)
    %histogram(alpha_all)
    
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

    scatter_twomoons_classify(data, u_avg(:, num_iterations), label_data)
    fname = 'print_runs/final_scatter.png';
    print('-r144','-dpng',fname);
    
    % plotVotes(Z, [371 280], 268:435)
    
    final_tau = tau_avg(num_iterations);
    final_alpha = alpha_avg(num_iterations);
    
    % correct_percent = count_correct(u_avg(:, num_iterations), set_neg, set_pos, [zeros(267,1) - 1; zeros(168,1) + 1]);
    correct_percent = count_correct(u_avg(:, num_iterations), set_neg, set_pos, ...
        [zeros(N/2,1) + 1; zeros(N/2, 1) - 1]);
    
    run_description = 'MCMC self-tuning Laplacian, xi tau alpha parameterization.\n';
    
    print_info_file(run_description, num_iterations, burn_in, p, q, l, B, ...
        gamma, final_tau, final_alpha, set_neg, set_pos, correct_percent, tau_epsilon, ...
        alpha_epsilon, elapsed_time, init_tau, init_alpha);

end

function p = count_correct(final_avg, set_neg, set_pos, correct_labels)
p = 0;
remainder = length(correct_labels) - length(set_neg) - length(set_pos);

for i=1:length(correct_labels)
    if ~ismember(i, set_neg) && ~ismember(i, set_pos)
        p = p + 1 - abs(sign(final_avg(i)) - correct_labels(i))/2;
    end
end
p = p / remainder;

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

function print_info_file(run_description, num_iterations, burn_in, p, q, l, B, ...
    gamma, final_tau, final_alpha, set_neg, set_pos, correct_percent, tau_epsilon,...
    alpha_epsilon, elapsed_time, init_tau, init_alpha)
fileID = fopen('print_runs/run_info.txt','w');
fprintf(fileID, 'This is an autogenerated text file with run info and summary!\n');
fprintf(fileID, ['Description: ' run_description]);
fprintf(fileID, 'Iterations = %d, Burn in period = %d\n', num_iterations, burn_in);
fprintf(fileID, 'Parameters of weight function: p = %d, q = %d, l = %d\n', p, q, l);
fprintf(fileID, 'Beta = %f, Gamma = %d\n', B, gamma);
fprintf(fileID, 'Tau epsilon = %f, Alpha epsilon = %f\n', tau_epsilon, alpha_epsilon);
fprintf(fileID, 'Initial tau = %f, Initial alpha = %f\n', init_tau, init_alpha);
fprintf(fileID, 'Average tau = %f, Average alpha = %f\n', final_tau, final_alpha);
fprintf(fileID, 'Label Data:\n');
fprintf(fileID, '+1: %d\n', set_pos);
fprintf(fileID, '-1: %d\n', set_neg);
fprintf(fileID, 'Percent correctly classified: %f\n', correct_percent);
fprintf(fileID, 'Time elapsed: %.2f s\n', elapsed_time);

end

function plotAvg(f, num_iterations, plot_me, k)
    figure(k)
    num_plots = length(plot_me);
    for i = 1:num_plots
        subplot(2, num_plots/2, i);
        plot(1:num_iterations, f(plot_me(i), :))
        title(sprintf('Senator %d', plot_me(i)))
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Running Averages','HorizontalAlignment','center','VerticalAlignment', 'top');
end

function plotVar(var, num_iterations, plot_me, k)
    figure(k)
    num_plots = length(plot_me);
    for i = 1:num_plots
        subplot(2, 4, i);
        plot(1:num_iterations, var(plot_me(i), :))
        title(sprintf('Senator %d', plot_me(i)))
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Running Variances','HorizontalAlignment','center','VerticalAlignment', 'top');
    
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
