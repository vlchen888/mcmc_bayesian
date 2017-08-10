fileID = fopen('MNIST_tests/record.csv','a');
fmt = "Trial #,Method,Accuracy\n";
fprintf(fileID,fmt);

for i = 1:5
    rng(i)
    correct_p = test_mcmc_multiclass_t_a_M_same(.01);
    fprintf(fileID,"%d,Hierarchical,",i);
    for j=1:length(correct_p)
        fprintf(fileID,"%.4f",correct_p(j));
        if j~=length(correct_p)
            fprintf(fileID,",");
        end
    end
    fprintf(fileID,"\n");

    rng(i)
    correct_p = test_mcmc_multiclass_TEST(.01);
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