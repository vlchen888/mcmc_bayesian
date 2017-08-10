rng(1123);test_mcmc_multiclass_t_a_M_same(.03);
params('laplacian') = "self tuning";
digs = [1,3,4,5,9];

Sample number: 100000, Time elapsed: 1032.52
Classification accuracy with S(E(u)): 0.9453
Classification accuracy with S(E(S(u))): 0.9448
    xi accept acceptance probability: 0.4532
    xi accept acceptance probability: 0.5639
    xi accept acceptance probability: 0.1227
    xi accept acceptance probability: 0.3590
    xi accept acceptance probability: 0.1282
    tau accept acceptance probability: 0.2170
    alpha accept acceptance probability: 0.5989
    M accept acceptance probability: 0.0000
tau running average: 0.7208
alpha running average: 53.4597

params('num_iterations') = 100001;
burn_in = 5000;

params('gamma') = 0.0001;
params('B') = 0.1;

params('init_tau') = 1.5;
params('tau_epsilon') = 0.08;
params('tau_min') = 0.01;
params('tau_max') = 20;


params('init_alpha') = 35;
params('alpha_epsilon') = 2;
params('alpha_min') = 0.1;
params('alpha_max') = 90;

params('init_M') = 50;
params('M_max_jump') = 0;
params('M_min') = 1;
params('M_max') = 50;