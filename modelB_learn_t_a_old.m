% mcmc_cluster();


function [tau_all, alpha_all, std] = mcmc_cluster()
    start_time = cputime;
    num_iterations = 1000;
    burn_in = 1;
    
    load('data3.mat')

    p = 2;
    q = 2;
    l = 1;

    gamma       = 1;
    B           = 0.4;
    Z           = X;
    set_neg     = [25 26 27 28];
    set_pos     = [280 281 282 283];
    
    init_tau        = 20;
    init_alpha      = 1;
    
    min_tau         = 0;
    max_tau         = 60;
    
    min_alpha       = 0;
    max_alpha       = 10;
    
    alpha_epsilon   = 0.1; % Jump alpha
    tau_epsilon     = 1; % Jump tau
     
    [L, ~, ~] = computeLaplacian(Z, p, q, l);
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    [num_senators, ~] = size(Z);
    
    label_data = init(num_senators, set_neg, set_pos);
    
    U = zeros(num_senators, num_iterations);
    
    %%%%% Indicates initialization from Fiedler Vector %%%%%
    U(2, 1) = 1;
            
    tau_all = zeros(1, num_iterations);
    tau_all(1) = init_tau;
    
    alpha_all = zeros(1, num_iterations);
    alpha_all(1) = init_alpha;
    
    %%%%% Acceptance probabilities %%%%%
    u_accept        = zeros(1, num_iterations);
    tau_accept      = zeros(1, num_iterations);
    alpha_accept    = zeros(1, num_iterations);
    
    u_accept_avg        = zeros(1, num_iterations);
    tau_accept_avg      = zeros(1, num_iterations);
    alpha_accept_avg    = zeros(1, num_iterations);
    
    %%%%% Take averages of tau, alpha over time as well %%%%%
    tau_avg = zeros(1, num_iterations);
    alpha_avg = zeros(1, num_iterations);
    
    %%%%%%% Do running average and variance as well %%%%%%
    % f(j, n) represents the average of the sign of the first n clusterings
    f = zeros(num_senators, num_iterations);
    var = zeros(num_senators, num_iterations);
    
    %%%%% Store standard basis vectors %%%%%
    std = zeros(num_senators, num_iterations);
    std(:, 1) = convert_std_basis(U(:, 1), phi);
    
    % Initialize f
    % for j=1:num_senators
    %    f(j, 1)=sign(u_0(j));
    % end
    
    for i=1:num_iterations-1
        
        %%%%% Propose new state for U %%%%%
        tau = tau_all(i);
        alpha = alpha_all(i);
        x = compute_rand_in_eigenbasis(lambda, tau, alpha);
        V_eigenbasis = (1-B^2)^0.5*U(:,i)+B*x;
        V = convert_std_basis(V_eigenbasis, phi);
        
        %%%%% Compute Transition Probability U -> V %%%%%
        transition_uv = min(1, likelihood(gamma, label_data, V)...
            /likelihood(gamma, label_data, std(:, i)));
        
        %%%% Do transition %%%%
        if rand(1) < transition_uv 
            U(:, i+1) = V_eigenbasis;
            std(:, i+1) = V;
            u_accept(i+1) = 1;
        else
            U(:, i+1) = U(:, i);
            std(:, i+1) = std(:, i);
            u_accept(i+1) = 0;
        end
                
        %%%%% Propose a new tau %%%%%
        new_tau = tau_all(i) + tau_epsilon * normrnd(0, 1);
        if new_tau < min_tau
            new_tau = min_tau;
        elseif new_tau > max_tau
            new_tau = max_tau;
        end
        
        log_tau = compute_log_prior(U(:, i+1), lambda, new_tau, alpha_all(i)) ...
            - compute_log_prior(U(:, i+1), lambda, tau_all(i), alpha_all(i));
        transition_tau = exp(log_tau);
        if rand(1) < transition_tau
            tau_all(i+1) = new_tau;
            tau_accept(i+1) = 1;
        else
            tau_all(i+1) = tau_all(i);
            tau_accept(i+1) = 0;
        end
        
        %%%%% Propose a new alpha %%%%%
        new_alpha = alpha_all(i) + alpha_epsilon * normrnd(0, 1);
        
        if new_alpha < min_alpha
            new_alpha = min_alpha;
        elseif new_alpha > max_alpha
            new_alpha = max_alpha;
        end
                
        log_alpha = compute_log_prior(U(:, i+1), lambda, tau_all(i+1), new_alpha) ...
            - compute_log_prior(U(:, i+1), lambda, tau_all(i+1), alpha_all(i));
        transition_alpha = exp(log_alpha);
        if rand(1) < transition_alpha
            alpha_all(i+1) = new_alpha;
            alpha_accept(i+1) = 1;
        else
            alpha_all(i+1) = alpha_all(i);
            alpha_accept(i+1) = 0;
        end
        
        %%%%% Update average and variance %%%%%
        if i >= burn_in
            f(:, i+1) = ((i-burn_in)*f(:, i) + sign(std(:, i+1)))/(i-burn_in+1);
            var(:, i+1) = 1-f(:, i+1).^2;
            
            tau_avg(i+1) = ((i-burn_in)*tau_avg(i) + tau_all(i+1))/(i-burn_in+1);
            alpha_avg(i+1) = ((i-burn_in)*alpha_avg(i) + alpha_all(i+1))/(i-burn_in+1);
            
            u_accept_avg(i+1) = ((i-burn_in)*u_accept_avg(i) + u_accept(i+1))/(i-burn_in+1);
            tau_accept_avg(i+1) = ((i-burn_in)*tau_accept_avg(i) + tau_accept(i+1))/(i-burn_in+1);
            alpha_accept_avg(i+1) = ((i-burn_in)*alpha_accept_avg(i) + alpha_accept(i+1))/(i-burn_in+1);
        end
        
        %%%%%% CODE FOR AVG MOVIE %%%%
        %if mod(i, 100) == 1   
        %    plotBar(f(:, i))
        %    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        %    str = sprintf('\\bf Step %d', i);
        %    text(0.5, 1, str,'HorizontalAlignment','center','VerticalAlignment', 'top');
        %    fname = sprintf('figs/step_%i.png',i);
        %    print('-r144','-dpng',fname);
        %end
    end
    
    plot_me = [(1:10)'; (268:277)'];

    %%%%% PRINT RUNNING AVERAGE OF U %%%%%
    %plotAvg(f, num_iterations, plot_me, 1)
    %fname = 'print_runs/running_avg.png';
    %print('-r144','-dpng',fname);
    
    %%%%% PRINT TRACES OF TAU AND ALPHA %%%%%
    clf
    plot(1:num_iterations, tau_all);
    ylabel('\tau trace');
    fname = 'print_runs/tau_all.png';
    print('-r144','-dpng',fname);
    
    clf
    plot(1:num_iterations, alpha_all);
    ylabel('\alpha trace');
    fname = 'print_runs/alpha_all.png';
    print('-r144','-dpng',fname);
    
    %%%%%%% PRINT RUNNING VARIANCES %%%%%
    %plotVar(var, num_plots, num_iterations, plot_me)
    
    %%%%% PRINT CLUSTERING %%%%%
    %plotCluster(f, num_iterations, num_senators)
    
    %%%%%%% FINAL VARIANCE %%%%%%%
    %plotBar(var(:, num_iterations))
    
    %%%%%%% FINAL AVG %%%%%%%
    plotBar(f(:, num_iterations), 2)
    fname = 'print_runs/final_avg.png';
    print('-r144','-dpng',fname);
    
    %%%%%%% Plot Tau? %%%%%%
    %figure(3)
    %histogram(tau_all)
    
    %figure(4)
    %histogram(alpha_all)
    
    figure(5)
    plot(1:num_iterations, u_accept_avg)
    ylabel('u acceptance probability');
    fname = 'print_runs/u_acceptance_probability.png';
    print('-r144','-dpng',fname);
    
    figure(6)
    plot(1:num_iterations, tau_accept_avg)
    ylabel('\tau acceptance probability');
    fname = 'print_runs/tau_acceptance_probability.png';
    print('-r144','-dpng',fname);
    
    figure(7)
    plot(1:num_iterations, alpha_accept_avg)
    ylabel('\alpha acceptance probability');
    fname = 'print_runs/alpha_acceptance_probability.png';
    print('-r144','-dpng',fname);

    
    %plotVotes(Z, [371 280], 268:435)
    
    final_tau = tau_avg(num_iterations);
    final_alpha = alpha_avg(num_iterations);
    correct_percent = count_votes_correct(f(:, num_iterations), set_neg, set_pos, [zeros(267,1) - 1; zeros(168,1) + 1]);

    run_description = 'MCMC with unnormalized Laplacian prior, learning alpha and tau, and nonzero gamma.\n';
    
    elapsed_time = cputime - start_time;
    
    print_info_file(run_description, num_iterations, burn_in, p, q, l, B, ...
        gamma, final_tau, final_alpha, set_neg, set_pos, correct_percent, tau_epsilon, ...
        alpha_epsilon, elapsed_time, init_tau, init_alpha);

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
fprintf(fileID, 'Final tau = %f, Final alpha = %f\n', final_tau, final_alpha);
fprintf(fileID, 'Label Data:\n');
fprintf(fileID, '+1: %d\n', set_pos);
fprintf(fileID, '-1: %d\n', set_neg);
fprintf(fileID, 'Percent of senators correctly classified: %f\n', correct_percent);
fprintf(fileID, 'Time elapsed: %.2f s\n', elapsed_time);

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

function log_prior = compute_log_prior(x, lambda, tau, alpha)
    log_prior = -0.5*compute_log_det_c(lambda, tau, alpha) ...
        -0.5*compute_inner_prod(x, lambda, tau, alpha);
end

function prior = computeRatioPrior(x, lambda, tau, alpha)
    detc = computeDetC(lambda, tau, alpha);
    inprod = compute_inner_prod(x, lambda, tau, alpha);
    prior = 1/sqrt(detc) ...
        * exp(-0.5*inprod);
end

function innerProd = compute_inner_prod(u, lambda, tau, alpha)
    innerProd = 0;
    for j = 1:length(lambda)
        innerProd = innerProd + u(j)^2 * (lambda(j) + tau^2)^alpha;
    end
end

function log_det_c = compute_log_det_c(lambda, tau, alpha)
    sum = 0;
    for j = 1:length(lambda)
        sum = sum + log(lambda(j) + tau^2);
    end
    log_det_c = -alpha * sum;
end

function detC = computeDetC(lambda, tau, alpha)
    prod = 1;
    for i = 2:length(lambda)
        prod = prod * (lambda(i)+ tau^2)^-alpha;
    end
    detC = prod;
end

function p = count_votes_correct(final_avg, set_neg, set_pos, correct_labels)
p = 0;
remainder = length(correct_labels) - length(set_neg) - length(set_pos);

for i=1:length(correct_labels)
    if ~ismember(i, set_neg) && ~ismember(i, set_pos)
        p = p + 1 - abs(sign(final_avg(i)) - correct_labels(i))/2;
    end
end
p = p / remainder;
end

function l = likelihood(gamma, label_data, u)
    sum = 0;
    for i=1:length(label_data)
        %%%% Check if i \in Z' %%%%
        if label_data(i) ~= 0
            sum = sum + abs(sign(u(i)) - label_data(i))^2;
        end
    end
    l = exp(-sum/(2*gamma^2));
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

function x = compute_rand_in_eigenbasis(lambda, tau, alpha)
    x = zeros(length(lambda), 1);
    for j=1:length(lambda)
       zeta = normrnd(0, 1);
       x(j) = eigvalweight(lambda(j), tau, alpha) * zeta;
    end
end

function u = convert_std_basis(x, phi)
    u = 0;
    for j=1:length(phi)
        u = u + x(j)* phi(:,j);
    end
end

function w = eigvalweight(lambda, tau, alpha)
    w = (lambda + tau^2)^(-alpha/2);
end

function plotEigenvectors(L)
    [V, ~] = eig(L);
    
    % Plot eigenvectors and determine clustering
    num_plots = 1;
    figure(4)
%     for i = 1:num_plots
%         subplot(2, 2, i);
%         plot(1:length(V(:, i)), V(:, i))
%     end
    plot(1:length(V(:,2)), V(:, 2))

end

function [L, L_sym, L_rw] = computeLaplacian(X, p, q, l)
    %compute Laplacian matrices
    [n, ~] = size(X);
    
    W = zeros(n);
    D = zeros(n);
    for i=1:n
        sum = 0;
        for j=1:n
            x = X(i,:);
            y = X(j,:);
            W(i, j) = weight(x, y, l, p, q);
            sum = sum + W(i, j);
        end
        D(i, i) = sum;
    end
    L = D - W; 
    L_sym = sqrtm(inv(D)) * L * sqrtm(inv(D));
    L_rw = inv(D) * L;
   %plotEigenvectors(L)
end

function d = dist(x, y, p, q)
%Lp of x - y raised to the qth power
    d = norm(x - y, p)^q;
end

function w = weight(x, y, l, p, q)
%weight function with given metric and length-scale parameter
    w = exp(-dist(x, y, p, q)/ 2 / l^2);
end