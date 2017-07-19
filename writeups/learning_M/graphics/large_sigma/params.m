percent_fielity = .03;
sigma = 0.2;

params('data_set') = string('moons');
params('laplacian') = string('self tuning');

params('num_iterations') = 100000;
burn_in = 5000;
params('burn_in') = burn_in;

params('gamma') = 0.1;
params('B') = 0.05;
params('init_tau') = 1;
params('init_alpha') = 35;
params('init_M') = 30;

params('min_tau')       = 0.1;
params('max_tau')       = 60;

params('min_alpha')     = 0.1;
params('max_alpha')     = 60;

params('min_M')       = 1;
params('max_M')       = 200;

params('alpha_epsilon') = 0;
params('tau_epsilon')   = 0;

params('max_M_jump') = 7;

p = 0.9021;