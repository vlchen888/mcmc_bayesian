rng(1111)

mnist 4 9
laplacian self tuning
params('num_iterations') = 100001;
burn_in = 5000;
params('burn_in') = burn_in;
params('movie') = true;

params('gamma') = 0.0001;
params('B') = 0.1;

params('init_tau')      = 10;
params('tau_epsilon')   = 0.1;
params('min_tau')       = 0.01;
params('max_tau')       = 20;

params('init_alpha')    = 35;
params('alpha_epsilon') = 0.3;
params('min_alpha')     = 0.1;
params('max_alpha')     = 60;

params('init_M')        = 50;
params('max_M_jump')    = 0;
params('min_M')         = 1;
params('max_M')         = 50;