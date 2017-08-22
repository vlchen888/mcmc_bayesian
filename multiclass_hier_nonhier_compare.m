function multiclass_hier_nonhier_compare()
%Comparing multiclass hierarchical and nonhierarchical algorithms
%   1% fidelity, gamma = 0.0001, MNIST [3, 4, 5, 9] tests, 30000 iterations
    fileID = fopen('MNIST_tests/record.csv','a');
    fprintf(fileID,"Comparing multiclass algorithms, .01 fidelity, gamma = 0.0001, MNIST3459, 30000 iterations\n");
    fprintf(fileID,"Trial #,Method,Sample,Accuracy\n");
    for i=8:15
        tic;
        rng(i^2);
        cont = test_mcmc_multiclass_t_a_M_same(.01);
        accuracy_map = cont('accuracy_map');
        
        fprintf(fileID,"%d,Hierarchical,\n",i);
        for sample = cell2mat(accuracy_map.keys)
            fprintf(fileID,",,%d,%.4f,\n",sample, accuracy_map(sample));
        end
        fprintf('Trial %d, hierarchical. Elapsed time: %.4f\n',i,toc);
        
        tic;
        rng(i^2);
        cont = test_mcmc_multiclass_simple(.01);
        accuracy_map = cont('accuracy_map');
        
        fprintf(fileID,"%d,Nonhierarchical,\n",i);
        for sample = cell2mat(accuracy_map.keys)
            fprintf(fileID,",,%d,%.4f,\n",sample, accuracy_map(sample));
        end
        fprintf('Trial %d, nonhierarchical. Elapsed time: %.4f\n',i,toc);
    end
    fclose(fileID);
end