function multiclass_hier_nonhier_compare()
%Comparing multiclass hierarchical and nonhierarchical algorithms
%   1% fidelity, gamma = 0.0001, MNIST [3, 4, 5, 9] tests, 30000 iterations
    fileID = fopen('MNIST_tests/record.csv','a');
    fmt = "Trial #,Method,Accuracy\n";
    fprintf(fileID,fmt);
    for i=1:15
        rng(i^2);
        cont = test_mcmc_multiclass_t_a_M_same(.01);
        correct_p = cont('correct_p');
        
        fprintf(fileID,"%d,Hierarchical,",i);
        for j=1:length(correct_p)
            fprintf(fileID,"%.4f",correct_p(j));
            if j~=length(correct_p)
                fprintf(fileID,",");
            end
        end
        fprintf(fileID,"\n");
        
        rng(i^2);
        cont = test_mcmc_multiclass_simple(.01);
        correct_p = cont('correct_p');
        
        fprintf(fileID,"%d,Nonhierarchical,",i);
        for j=1:length(correct_p)
            fprintf(fileID,"%.4f",correct_p(j));
            if j~=length(correct_p)
                fprintf(fileID,",");
            end
        end
        fprintf(fileID,"\n");
    end
    fclose(fileID);
end

