rng(5)

percent_fielity = .0115;

params('data_set') = string('voting');
params('laplacian') = string('un');

params('num_iterations') = 100000;
burn_in = 5000;
params('gamma') = 0.1;
params('B') = 0.4;
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

p = 0.8744;

-1:
90
97
213

+1:
379
400