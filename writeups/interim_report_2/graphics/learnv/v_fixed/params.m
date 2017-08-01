rng(8)

percent_fielity = .01;
sigma = 0.2;

params('num_iterations') = 100000;
burn_in = 5000;

params('gamma') = 0.1;
params('B') = 0.1;
params('tau') = 2;
params('alpha') = 35;
params('init_M') = 50;
params('min_M')       = 1;
params('max_M')       = 70;
params('max_M_jump') = 20;

p = 0.8652;
