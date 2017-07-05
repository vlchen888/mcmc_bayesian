function p = test_mcmc_hier_centered()
params = containers.Map;

load('data3.mat');
data = X;
params('data') = data;
params('num_iterations') = 10000;
burn_in = 1000;
params('p') = 2;
params('q') = 2;
params('l') = 1.25;

params('gamma')         = 0.0001;
params('B')             = 0.6;
params('init_tau')      = 20;
params('init_alpha')    = 20;

params('min_tau')       = 0.1;
params('max_tau')       = 60;

params('min_alpha')     = 0.1;
params('max_alpha')     = 60;

params('alpha_epsilon') = 0.1;
params('tau_epsilon')   = 0.1;

params('data_set') = string('voting');
params('laplacian') = string('un');

if params('data_set') == string('voting')
    load('data3.mat')
    data = X;
    params('data') = data;
    set_neg     = 20:22;
    set_pos     = 280:281;
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

[tau_all, alpha_all, std, u_accept, tau_accept, alpha_accept] = mcmc_learn_t_a(params);
u_avg = mean(sign(std(:, burn_in:end)), 2); %avg the rows
figure(1)
clf
plotBar(u_avg);
p = count_correct(u_avg, label_data, [zeros(267,1) - 1; zeros(168,1) + 1]);
u_avg_accept = mean(u_accept(burn_in:end))
tau_avg_accept = mean(tau_accept(burn_in:end))
alpha_avg_accept = mean(alpha_accept(burn_in:end))
end

function u = init(num_data, set_neg, set_pos)
    u = zeros(num_data, 1);
    
    % label some of the data
    for i=1:length(set_pos)
        u(set_pos(i))=1;
    end
    for i=1:length(set_neg)
        u(set_neg(i))=-1;
    end
end
