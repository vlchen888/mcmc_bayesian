rng(1123);test_mcmc_multiclass_TEST(.03);
params('laplacian') = "self tuning";
digs = [1,3,4,5,9];

Sample number: 100000, Time elapsed: 864.35
Classification accuracy with S(E(u)): 0.9502
Classification accuracy with S(E(S(u))): 0.9502
    xi accept acceptance probability: 0.4409
    xi accept acceptance probability: 0.6129
    xi accept acceptance probability: 0.3938
    xi accept acceptance probability: 0.4328
    xi accept acceptance probability: 0.3681
    tau accept acceptance probability: 0.0000
    alpha accept acceptance probability: 0.0000
    M accept acceptance probability: 0.0000
tau running average: 0.0000
alpha running average: 1.0000

params('num_iterations') = 100001;
burn_in = 5000;

params('gamma') = 0.0001;
params('B') = 0.1;

params('init_tau') = 0;
params('tau_epsilon') = 0;
params('tau_min') = 0.01;
params('tau_max') = 20;


params('init_alpha') = 1;
params('alpha_epsilon') = 0;
params('alpha_min') = 0.1;
params('alpha_max') = 90;

params('init_M') = 50;
params('M_max_jump') = 0;
params('M_min') = 1;
params('M_max') = 50;