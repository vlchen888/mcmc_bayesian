rng(5)

percent_fielity = .0115;

params('data_set') = string('voting');
params('laplacian') = string('un');

params('num_iterations') = 100000;
burn_in = 5000;
params('gamma') = 0.0001;
params('B') = 0.4;
params('init_tau') = 2;
params('init_alpha') = 35;
M = 435;
p = 0.8767;
