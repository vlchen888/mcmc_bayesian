rng(3)

percent_fielity = .03;
sigma = 0.2;

params('data_set') = string('moons');
params('laplacian') = string('self tuning');

params('num_iterations') = 100000;
burn_in = 5000;
params('burn_in') = burn_in;
params('gamma') = 0.1;
params('B') = 0.1;
params('init_tau') = 1;
params('init_alpha') = 35;

p = 0.8567;
