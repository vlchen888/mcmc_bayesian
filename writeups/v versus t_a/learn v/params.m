rng(1111)

percent_fielity = .01;
sigma = 0.2;
params('data_set') = "moons";
params('laplacian') = "self tuning";

params('num_iterations') = 100000;
burn_in = 5000;

params('gamma') = 0.1;
params('B') = 0.1;
params('a') = 0.8;

params('epsilon') = .01;
params('tau') = 3 or 5;
params('alpha') = 35;

params('init_M') = 50;
params('max_M_jump') = 0;
params('min_M') = 1;
params('max_M') = 50;


tau = 3:
Sample number: 100000, Elapsed time: 107.0941
Classification accuracy S(E(S(u))): 0.861616
Acceptance rates: 
    xi: 0.6386
    v: 0.6378
    M: 0.0000


tau = 5:
Sample number: 100000, Elapsed time: 111.4390
Classification accuracy S(E(S(u))): 0.804545
Acceptance rates: 
    xi: 0.5911
    v: 0.6150
    M: 0.0000