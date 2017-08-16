rng(1111)

mnist 4 9
laplacian self tuning
params('num_iterations') = 100001;
burn_in = 5000;
params('burn_in') = burn_in;

params('gamma') = 0.0001;
params('B') = 0.1;
params('a') = 0.8;

%params('epsilon') = 1e-13;

params('epsilon') = .01;
params('tau') = 1;
params('alpha') = 35;

params('init_M') = 50;
params('max_M_jump') = 0;
params('min_M') = 1;
params('max_M') = 50;