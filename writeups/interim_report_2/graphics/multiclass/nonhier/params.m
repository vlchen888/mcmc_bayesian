rng(1234)
digs = [1, 4, 9];
params('laplacian') = "self tuning";
params('gamma') = 0.1;
params('B') = 0.1;
params('tau') = 2;
params('alpha') = 35;
params('num_iterations') = 100001;
burn_in = 5000;

Sample number: 100000, Time elapsed: 244.21
Classification accuracy with S(E(u)): 0.8964
Classification accuracy with S(E(S(u))): 0.8867
    xi accept acceptance probability: 0.6449
    xi accept acceptance probability: 0.4492
    xi accept acceptance probability: 0.4474