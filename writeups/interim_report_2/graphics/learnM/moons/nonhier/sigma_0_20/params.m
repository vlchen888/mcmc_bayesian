rng(2)

percent_fielity = .01;
sigma = 0.2;

params('data_set') = string('moons');
params('laplacian') = string('self tuning');

params('num_iterations') = 100000;
burn_in = 5000;
params('burn_in') = burn_in;
params('gamma') = 0.1;
params('B') = 0.4;
params('init_tau') = 2;
params('init_alpha') = 35;

M = 100;

p = 0.9056;
