
% data
%
sigma = 0.03;
gamma = 0.1;
num_data = 2000;
data = moondata(1,100,num_data,sigma);

J1 = 10;
J2 = 10;
J = J1+J2;

J = round(0.03*num_data);
J1 = randsample(J,1);
J2 = J-J1;

%jump = 268;
jump = num_data/2;
truth = [-ones(jump-1,1);ones(num_data-jump+1,1)];

xMinus = randsample(jump-1,J1);
xPlus = randsample(num_data-jump,J2)+jump;

label_data = zeros(num_data,1);
label_data(xMinus) = -1;
label_data(xPlus) = +1;

L = compute_laplacian_selftuning(data);

[phi, lambda] = eig(L);
lambda = diag(lambda);

M = 50;
phi = phi(:,1:M);
lambda = lambda(1:M);
%



% other parameters

num_iterations = 100000;
p = 2;
q = 2;
l = 1;
%gamma = 0.01;
B = 0.4;
init_tau = 1;
init_alpha = 35;
min_tau = -50;
max_tau = 50;
min_alpha = 1;
max_alpha = 50;
alpha_epsilon = 0.5;
tau_epsilon = 0.1;

tic

U = zeros(M, num_iterations);

%%%%% Indicates initialization from Fiedler Vector %%%%%
U(1, 1) = (lambda(1) + init_tau^2)^(-init_alpha/2);
U(2, 1) = -(lambda(2) + init_tau^2)^(-init_alpha/2);

tau_all = zeros(1, num_iterations);
tau_all(1) = init_tau;

alpha_all = zeros(1, num_iterations);
alpha_all(1) = init_alpha;

%%%%% Acceptance probabilities %%%%%
u_accept        = zeros(1, num_iterations);
tau_accept      = zeros(1, num_iterations);
alpha_accept    = zeros(1, num_iterations);

%%%%% Store standard basis vectors %%%%%
std = zeros(num_data, num_iterations);
std(:, 1) = convert_std_basis(U(:, 1), phi);

for i=1:num_iterations-1

    %%%%% Propose new state for U %%%%%
    tau = tau_all(i);
    alpha = alpha_all(i);
    
    %tau = init_tau;
    %alpha = init_alpha;
    
    x = compute_rand_in_eigenbasis(lambda, tau, alpha);
    V_eigenbasis = (1-B^2)^0.5*U(:,i)+B*x;
    V = convert_std_basis(V_eigenbasis, phi);

    %%%%% Compute Transition Probability U -> V %%%%%
    log_uv = compute_log_likelihood(gamma, label_data, V) - ...
        compute_log_likelihood(gamma, label_data, std(:, i));
    transition_uv = exp(log_uv);

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
    
    %

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
    %
    
    if mod(i,2500)==0
        subplot(131);
        plot(sign(phi)*U(:,i));
        %axis([0,num_data,-2,2]);
        
        subplot(132);
        scatter(data(phi*U(:,i)>0,1),data(phi*U(:,i)>0,2));
        hold on;
        scatter(data(phi*U(:,i)<0,1),data(phi*U(:,i)<0,2));
        hold off;
        
        subplot(133);
        UMean = mean(U(:,1:i)')';
        scatter(data(phi*UMean>0,1),data(phi*UMean>0,2));
        hold on;
        scatter(data(phi*UMean<0,1),data(phi*UMean<0,2));
        hold off;
        pause(0.01);
        
        err = sign(phi*UMean)-truth;
        err([xMinus;xPlus]) = [];
        
        fprintf('Sample Number: %i\n',i);
        fprintf('\tAcceptance rates: U:%f, t:%f, a:%f\n',sum(u_accept)/i,sum(tau_accept)/i,sum(alpha_accept)/i);
        fprintf('\tClassification Accuracy: %f\n',(num_data-J-length(find(err~=0)))/(num_data-J))
    end
    

end

toc

function log_prior = compute_log_prior(u, lambda, tau, alpha)
    log_prior = -0.5*compute_log_det_c(lambda, tau, alpha) ...
        - 0.5*compute_inner_prod(u, lambda, tau, alpha);
end

function innerProd = compute_inner_prod(u, lambda, tau, alpha)
    innerProd = (u.^2)' * (lambda + tau^2).^alpha;
end

function log_det_c = compute_log_det_c(lambda, tau, alpha)
    log_det_c = -alpha * sum(log(lambda + tau^2));
end

function l = compute_log_likelihood(gamma, label_data, u)
    sum = norm((sign(u)-label_data).*abs(label_data))^2;
    l = -sum/(2*gamma^2);
end

function x = compute_rand_in_eigenbasis(lambda, tau, alpha)
    x = eigvalweight(lambda, tau, alpha) .* normrnd(0, 1, length(lambda), 1);
end

function u = convert_std_basis(x, phi)
    u = phi*x;
end

function w = eigvalweight(lambda, tau, alpha)
    w = (lambda + tau^2).^(-alpha/2);
end