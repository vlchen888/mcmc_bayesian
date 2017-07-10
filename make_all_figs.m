record = csvread('csvread/fiedler.csv', 1); % Start reading from row 1
sigma_arr = [0.02 0.04 0.06 0.08 0.1];
fidelity_arr = [0.5 1 3];
%fidelity_arr = [0.005 0.01 0.03];
num_trials = 50;
figure(1)
clf
make_figs(record, sigma_arr, fidelity_arr, 'r')

fidelity_arr = [0.005 0.01 0.03];
record = csvread('csvread/mcmc_gamma.csv', 1); % Start reading from row 1
make_figs(record, sigma_arr, fidelity_arr, 'k')

record = csvread('csvread/centered.csv', 1); % Start reading from row 1
make_figs(record, sigma_arr, fidelity_arr, 'g')

record = csvread('csvread/noncentered.csv', 1); % Start reading from row 1
make_figs(record, sigma_arr, fidelity_arr, 'b')

legend('Fiedler','Nonhierarchical', 'Centered', 'Noncentered');