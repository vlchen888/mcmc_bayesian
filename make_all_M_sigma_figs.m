record = csvread('M_sigma_tests/record.csv', 1); %start from row 1
sigma_arr = [0.04 0.08 0.12 0.16 0.20];
num_trials = 50;
figure(1)
set(gcf, 'Position', [100, 300, 1200, 300])
clf
hold on
for metric = 4:4
    plot_me = [];
    plot_err_neg = [];
    plot_err_pos = [];
    for i = 1:length(sigma_arr)
        sigma = sigma_arr(i);
        TF = record(:,2) == sigma;
        trials = record(TF,:);
        med = median(trials(:,metric));
        L = med - prctile(trials(:,metric), 25);
        R = prctile(trials(:,metric), 75) - med;
        plot_me = [plot_me; med];
        plot_err_neg = [plot_err_neg; L];
        plot_err_pos = [plot_err_pos; R];
    end
    errorbar(sigma_arr,plot_me,plot_err_neg,plot_err_pos, 'LineWidth', 2)
end

xlabel('\sigma')
legend('Classification accuracy')
%legend('M_{mean}', 'M_{median}', 'M_{min}')