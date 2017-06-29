load('data3.mat');
mcmc_cluster();

% p = 2;
% q = 2;
% l = 1;
% numclusters = 2;
% computeLaplacian(X, p, q, l);


function u = init(num_senators)
    u = zeros(num_senators, 1);
    
    % label some of the data
    set_pos = [25 26 27 28];
    set_neg = [280 281 282 283];
    for i=1:length(set_pos)
        u(set_pos(i))=1;
    end
    for i=1:length(set_neg)
        u(set_neg(i))=-1;
    end
end

function c = check(u, label_data, num_senators)
    s = sign(label_data);
    
    for i = 1:num_senators
        if s(i) ~= 0 && s(i)*u(i)<0
            c = false;
            return
        end
    end
    c = true;
end

function mcmc_cluster()
    load('data3.mat')

    p = 2;
    q = 2;
    l = 1;
    a = 1;
    tau = 0;
    alpha = 1;
    gamma = 0.5;

    Z = X;
    [L, ~, ~] = computeLaplacian(Z, p, q, l);
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    Z = X;
    [num_senators, ~] = size(Z); % uses columns as senators
    
    label_data = init(num_senators);
    
    u_0 = zeros(num_senators, 1);
    for i=1:num_senators
        if label_data(i) ~= 0
            u_0(i) = label_data(i);
        else
            u_0(i) = phi(i, 2);
        end
    end
            
    num_iterations = 1000;
    
    U = zeros(num_senators, num_iterations);
    U(:, 1) = u_0;
    
    B = 0.4;
    %%%%%%% Do running average and variance as well %%%%%%
    % f(j, n) represents the average of the sign of the first n clusterings
    f = zeros(num_senators, num_iterations);
    var = zeros(num_senators, num_iterations);
    
    % Initialize f
    for j=1:num_senators
        f(j, 1)=sign(u_0(j));
    end
    
    for i=1:num_iterations
        if i > 1
            u_star = (1-B^2)^0.5*U(:,i)+B*compute_rand(lambda, phi, num_senators, a, tau, alpha);

            if check(u_star, label_data, num_senators)
                U(:, i) = u_star;
            else
                U(:, i) = U(:, i-1);
            end

            f(:, i) = ((i-1)*f(:, i-1) + sign(U(:, i)))/i;
            var(:, i) = 1-f(:, i).^2;
        end
        
        %%%%%% CODE FOR AVG MOVIE %%%%
        %plotBar(f(:, i))
        %ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        %str = sprintf('\\bf Step %d', i);
        %text(0.5, 1, str,'HorizontalAlignment','center','VerticalAlignment', 'top');
        %fname = sprintf('figs/step_%i.png',i);
        %print('-r144','-dpng',fname);
    end
    
    plot_me = [21 22 23 1 435 434 433 421];
    num_plots = length(plot_me);

    %%%%%%% PRINT RUNNING AVERAGE %%%%%
    plotAvg(f, num_plots, num_iterations, plot_me)
    
    %%%%%%% PRINT RUNNING VARIANCES %%%%%
    plotVar(var, num_plots, num_iterations, plot_me)
    
    %%%%% PRINT CLUSTERING %%%%%
    %plotCluster(f, num_iterations, num_senators)
    
    %%%%%%% FINAL VARIANCE %%%%%%%
    %plotBar(var(:, num_iterations))
    
    %%%%%%% FINAL AVG %%%%%%%
    plotBar(f(:, num_iterations))
    
    plotVotes(Z, [371 280], 268:435)
    

    
end

function plotVotes(Z, plot_diff, plot_avg)
    [~, num_votes] = size(Z);
    figure(6)
    clf
    %plot_legend = string(zeros(length(plot_me), 1));
    num_plots = length(plot_diff) + 1;
    
    for i = 1:length(plot_diff)
        subplot(num_plots, 1, i)
        plotted = bar(1:num_votes, Z(plot_diff(i),:));
        legend(plotted, string(plot_diff(i)))
    end
    
    avg = zeros(1, num_votes);
    for i = 1:length(plot_avg)
        sen = plot_avg(i);
        avg = avg + Z(sen, :);
    end
    avg = avg / length(plot_avg);
    
    subplot(num_plots, 1, num_plots)
    plotted = bar(1:num_votes, avg);
    legend(plotted, 'Average');
end

function plotAvg(f, num_plots, num_iterations, plot_me)
    figure(1)
    
    for i = 1:num_plots
        subplot(2, 4, i);
        plot(1:num_iterations, f(plot_me(i), :))
        title(sprintf('Senator %d', plot_me(i)))
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Running Averages','HorizontalAlignment','center','VerticalAlignment', 'top');
end

function plotVar(var, num_plots, num_iterations, plot_me)
    figure(2)
    for i = 1:num_plots
        subplot(2, 4, i);
        plot(1:num_iterations, var(plot_me(i), :))
        title(sprintf('Senator %d', plot_me(i)))
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Running Variances','HorizontalAlignment','center','VerticalAlignment', 'top');
    
end

function plotBar(avg)
    
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
    
    
    figure(5)
    clf
    hold on
    bar(posx, posy, 'b')
    bar(negx, negy, 'r')
    hold off
end

function plotCluster(f, num_iterations, num_senators)
    clusters = sign(f(:, num_iterations));
    
    
    negative_x = zeros(435, 1);
    negative_y = zeros(435, 1);
    num_neg = 0;
    
    positive_x = zeros(435, 1);
    positive_y = zeros(435, 1);
    num_pos = 0;
    
    zero_x = zeros(435, 1);
    zero_y = zeros(435, 1);
    num_zero = 0;
    
    for i=1:num_senators
        if clusters(i) < 0
            num_neg = num_neg + 1;
            negative_x(num_neg) =  i; 
            negative_y(num_neg) = -1; 
            
        elseif clusters(i) > 0
            num_pos = num_pos + 1;
            positive_x(num_pos) =  i; 
            positive_y(num_pos) = 1; 
        else
            num_zero = num_zero + 1;
            negative_x(num_zero) = i;
            negative_y(num_zero) = 0;
        end
    end
    figure(3)
    clf
    hold on
    xlim([-inf inf])
    ylim([-2 2])
    scatter(positive_x(1:num_pos), positive_y(1:num_pos), 18, 'blue');
    scatter(negative_x(1:num_neg), negative_y(1:num_neg), 18, 'red');
    scatter(zero_x(1:num_zero), zero_y(1:num_zero), 18, 'black');
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Clustering','HorizontalAlignment','center','VerticalAlignment', 'top');

    hold off
end

function u = compute_rand(lambda, phi, num_senators, a, tau, alpha)
    u = 0;
    for j=2:num_senators
        zeta = normrnd(0, 1);
        u = u + eigvalweight(lambda(j), a, tau, alpha) * zeta * phi(:,j);
    end
end

function w = eigvalweight(lambda, a, tau, alpha)
    w = a*(lambda + tau^2)^(-1/2*alpha);
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
   plotEigenvectors(L)
end

function d = dist(x, y, p, q)
%Lp of x - y raised to the qth power
    d = norm(x - y, p)^q;
end

function w = weight(x, y, l, p, q)
%weight function with given metric and length-scale parameter
    w = exp(-dist(x, y, p, q)/ 2 / l^2);
end