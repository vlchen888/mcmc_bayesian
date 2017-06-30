function p = compute_mcmc_gamma_accuracy()
params = containers.Map;

load('data3.mat');
data = X;
params('data') = data;
params('num_iterations') = 1000;
burn_in = 1;
params('p') = 2;
params('q') = 2;
params('l') = 1;
params('gamma') = 0.0001;
params('B') = 0.4;
params('init_tau') = 0.1;
params('init_alpha') = 1;
params('data_set') = string('voting');
params('laplacian') = string('un');

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

[u, u_accept] = mcmc_gamma(params);
u_avg = mean(sign(u(:, burn_in:end)), 2); %avg the rows
plotBar(u_avg);
p = count_correct(u_avg, label_data, [zeros(267,1) - 1; zeros(168,1) + 1]);

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
