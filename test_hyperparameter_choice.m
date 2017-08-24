function test_hyperparameter_choice()
%Testing (tau, M) vs (alpha, M)
%   Detailed explanation goes here
fileID = fopen('hyperparameter_tests/non_hier_11-20.csv','a');
fprintf(fileID,"Hyperparameter choice, .01 fidelity, gamma = 0.0001, MNIST3459\n");
fprintf(fileID,"Trial #,Method,Sample,Accuracy\n");
%descriptions = ["tau_alpha", "tau_M", "alpha_M", "tau_alpha_M"];
descriptions = ["simple"];
for i = 11:20
    for descrip_id = 1:length(descriptions)
        tic;
        rng(i^2);
        cont = test_choose_params(descriptions(descrip_id));
        accuracy_map = cont('accuracy_map');
        print_test(fileID, descrip_id-1, accuracy_map, i);
        fprintf('Trial %d, %s. Elapsed time: %.4f\n',i, descriptions(descrip_id), toc);
        fig_max = 4;
        if descriptions(descrip_id) == "simple"
           fig_max = 3; 
        end
        for fig = 1:fig_max % FOR SIMPLE MCMC ONLY
            figure(fig)
            fname = sprintf('hyperparameter_tests/figs/%s/trial_%i_fig_%i.png',descriptions(descrip_id),i,fig);
            print('-r144','-dpng',fname);
        end
    end
end
fclose(fileID);
end

function cont = test_choose_params(description)
    if description == "tau_alpha"
        cont = test_multiclass_t_a(.01);
    elseif description == "tau_M"
        cont = test_multiclass_t_M(.01);
    elseif description == "alpha_M"
        cont = test_multiclass_a_M(.01);
    elseif description == "tau_alpha_M"
        cont = test_multiclass_t_a_M(.01);
    elseif description == "simple"
        cont = test_mcmc_multiclass_simple(.01);
    end
end

function print_test(fileID, description, accuracy_map, i)
%fprintf(fileID,"%d,%s,,\n",i,description);
for sample = cell2mat(accuracy_map.keys)
        fprintf(fileID,"%d,%d,%d,%.4f\n",i,description, sample, accuracy_map(sample));
end
end
