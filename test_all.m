function test_all()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

sigma_arr = [0.02 0.04 0.06 0.08 0.1];
A = length(sigma_arr);

fidelity_arr = [0.005 0.01 0.03];
B = length(fidelity_arr);

num_trials = 50;

fileID = fopen('tests/record.csv','w');
fmt = 'Trial #,Sigma,Fidelity,Accuracy,tau,alpha\n';
fprintf(fileID,fmt);
fclose(fileID);
    for i = 1:num_trials
        fprintf('Trial %d started.\n', i)
        for a = 1:A
            for b = 1:B
                sigma = sigma_arr(a);
                percent_fidelity = fidelity_arr(b);

                rng(i)
                p = test_mcmc_gamma_accuracy(percent_fidelity, sigma);
                
                fileID = fopen('tests/record.csv','a');
                fmt = '%d,%f,%f,%f,%f,%f\n';
                fprintf(fileID,fmt, i, sigma, percent_fidelity, p,0,1);
                fclose(fileID);
                
                %test_mcmc_hier_centered(percent_fidelity, sigma);
                %test_mcmc_hier_noncentered(percent_fidelity, sigma);
            end
        end
    end
end

