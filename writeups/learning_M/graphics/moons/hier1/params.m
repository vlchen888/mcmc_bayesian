rng(3)

percent_fielity = .01;
sigma = 0.2;

params('data_set') = string('moons');
params('laplacian') = string('self tuning');

params('num_iterations') = 100000;
burn_in = 5000;
params('burn_in') = burn_in;

params('gamma') = 0.1;
params('B') = 0.1;
params('init_tau') = 2;
params('init_alpha') = 35;
params('init_M') = 50;

params('min_tau')       = 0.1;
params('max_tau')       = 60;

params('min_alpha')     = 0.1;
params('max_alpha')     = 60;

params('min_M')       = 1;
params('max_M')       = 70;

params('alpha_epsilon') = 0;
params('tau_epsilon')   = 0;

params('max_M_jump') = 20;

p = 0.8495;
