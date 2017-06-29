load('data3.mat');
mcmc_cluster();

% p = 2;
% q = 2;
% l = 1;
% numclusters = 2;
% computeLaplacian(X, p, q, l);


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

function mcmc_cluster()
    num_iterations = 5000;
    
    load('data3.mat')

    a = 1;
    tau = 0;
    alpha = 1;
    gamma = 1;
    B = 0.9;
    Z = X;
    set_neg = [25 26 27 28];
    set_pos = [280 281 282 283];
    %set_neg = 1:267;
    %set_pos = [];
    
    burn_in = 2000;
    
    L = compute_laplacian_selftuning(X);
    lambda = eig(L);
    [phi, ~] = eig(L);
    
    [num_senators, ~] = size(Z);
    
    label_data = init(num_senators, set_neg, set_pos);
    
    u_0 = zeros(num_senators, 1);
    for i=1:num_senators
        if label_data(i) ~= 0
            u_0(i) = label_data(i);
        else
            u_0(i) = phi(i, 2);
        end
    end
    
    U = zeros(num_senators, num_iterations);
    U(:, 1) = u_0;
    
    u_accept = zeros(1, num_iterations);
    u_accept_avg = zeros(1, num_iterations);

    %%%%%%% Do running average and variance as well %%%%%%
    % f(j, n) represents the average of the sign of the first n clusterings
    f = zeros(num_senators, num_iterations);
    var = zeros(num_senators, num_iterations);
    
    % Initialize f
    % for j=1:num_senators
    %    f(j, 1)=sign(u_0(j));
    % end
    
    for i=1:num_iterations-1
        u_star = (1-B^2)^0.5*U(:,i)+B*compute_rand(lambda, phi, num_senators, a, tau, alpha);

        %%%%% COMPUTE TRANSITION PROBABILITY %%%%%
        transition_p = min(1, likelihood(gamma, label_data, u_star)/likelihood(gamma, label_data, U(:,i)));

        if rand(1) < transition_p %%%% TRANSITION %%%%
            U(:, i+1) = u_star;
            u_accept(i+1) = 1;
        else
            U(:, i+1) = U(:, i);
            u_accept(i+1) = 0;
        end
        
        if i >= burn_in
            f(:, i+1) = ((i-burn_in)*f(:, i) + sign(U(:, i+1)))/(i-burn_in+1);
            var(:, i+1) = 1-f(:, i+1).^2;
            u_accept_avg(i+1) = ((i-burn_in)*u_accept_avg(i) + u_accept(i+1))/(i-burn_in+1);
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
    
    figure(1)
    plot(1:num_iterations, u_accept_avg)
    
    plot_me = [(1:10)'; (268:277)'];

    %%%%%%% PRINT RUNNING AVERAGE %%%%%
    plotAvg(f, num_iterations, plot_me, 2)
    
    %%%%%%% PRINT RUNNING VARIANCES %%%%%
    %plotVar(var, num_iterations, plot_me)
    
    %%%%% PRINT CLUSTERING %%%%%
    %plotCluster(f, num_iterations, num_senators)
    
    %%%%%%% FINAL VARIANCE %%%%%%%
    %plotBar(var(:, num_iterations))
    
    %%%%%%% FINAL AVG %%%%%%%
    plotBar(f(:, num_iterations), 3)
    
    %plotVotes(Z, [371 280], 268:435)
    
    count_votes_correct(f(:, num_iterations), set_neg, set_pos, [zeros(267,1) - 1; zeros(168,1) + 1])


    
end

function l = likelihood(gamma, label_data, u)
    sum = 0;
    for i=1:length(label_data)
        if label_data(i) ~= 0 %%%% TEST IN Z' %%%%
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

function L = compute_laplacian_selftuning(data)

lap = 'un';                 % Laplacian: 'un','sym' or 'rw'
dist_type = 'euclidean';    % Distance for weights
q = 2;                      % Power of distance for weights
ell = 1;                    % Length-scale for weights
K = 7;                     % Number of nearest neighbours

N = size(data,1);
W=zeros(N,N); 

[IDX,DX] = knnsearch(data,data,'K',K,'Distance',dist_type,'IncludeTies',true);

for i=1:N
    for j=1:length(IDX{i})
        l = IDX{i}(j);
        W(i,l) = exp(-DX{i}(j)^q/(2*ell^q*DX{i}(K)*DX{j}(K)));
        %W(i,l) = exp(-DX{i}(j)^q/(2*ell^2));
    end
end
% Take care of 0/0 cases and symmetrize, just in case
W(isnan(W)) = 1;
W = (W+W')/2;

D=diag(sum(W));
DD=diag(1./sqrt(sum(W)));
L = D-W;

if strcmp(lap,'sym')
    L = DD*L*DD;
elseif strcmp(lap,'rw')
    L = DD*DD*L;
end
L = (L+L')/2;
plotEigenvectors(L)
end


