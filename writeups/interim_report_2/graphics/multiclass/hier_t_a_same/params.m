rng(1234)
digs = [1, 4, 9];
params('laplacian') = "self tuning";
params('gamma') = 0.1;
params('B') = 0.1;
params('tau') = 2;
params('alpha') = 35;
params('num_iterations') = 100001;
burn_in = 5000;

Sample number: 100000, Time elapsed: 366.09
Classification accuracy with S(E(u)): 0.8862
Classification accuracy with S(E(S(u))): 0.9433
    xi accept acceptance probability: 0.8565
    xi accept acceptance probability: 0.2495
    xi accept acceptance probability: 0.2497
    tau accept acceptance probability: 0.3031
    alpha accept acceptance probability: 0.8908

tau running average: 0.6377
alpha running average: 37.3116