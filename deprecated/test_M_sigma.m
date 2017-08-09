function test_M_sigma()
%TESTS M AND SIGMA RELATION
    %fileID = fopen('M_sigma_tests/record_part2.csv','w');
    %fmt = 'Trial #,Sigma,Fidelity,Accuracy,M_mean,M_median,M_min\n';
    %fprintf(fileID,fmt);
    %fclose(fileID);

    sigma_arr = [0.04, 0.08, 0.12, 0.16, 0.2];
    fidelity = 0.01;
    for i = 33:50
        fprintf('Trial %d started.\n', i)
        tic;
        for s = 1:length(sigma_arr)
            curr_sigma = sigma_arr(s);
            rng(i)
            [p, M_mean, M_median, M_min] = test_mcmc_t_a_M(fidelity, curr_sigma);
            fileID = fopen('M_sigma_tests/record.csv','a');
            fmt = "%d,%.2f,%.2f,%.4f,%.4f,%d,%d\n";
            fprintf(fileID, fmt, i, curr_sigma, fidelity, p, M_mean, M_median, M_min);
            fclose(fileID);
        end
        toc
    end
end

