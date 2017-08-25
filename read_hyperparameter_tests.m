function read_hyperparameter_tests()
%Reading hyperparameter test results
methods = {'simple';'tau_alpha'; 'tau_M'; 'alpha_M'; 'tau_alpha_M'}; % types 0,1,2,3,4
record = csvread('hyperparameter_tests/all1-30.csv', 2);
samples = [5000 10000 15000 20000 25000];

for i = 1:length(methods)
    method = i - 1;
    set(gcf, 'Position', [0, 1000, 600, 400])
    plot_arr = zeros(length(samples),1);
    L_arr = zeros(length(samples),1);
    R_arr = zeros(length(samples),1);
    for sample = samples
        T = record(:,2) == method & record(:,3) == sample;
        these_records = record(T,4);
        plot_arr(sample/5000) = median(these_records);
        
        L_arr(sample/5000) = plot_arr(sample/5000) - prctile(these_records, 25);
        
        R_arr(sample/5000) = prctile(these_records, 75) - plot_arr(sample/5000);
    end
    plot(samples,plot_arr,'LineWidth',2);
    %errorbar(samples,plot_arr,L_arr,R_arr,'LineWidth',2);
    hold on;
end

title("Comparing different hyperparameter choices, 30 trials");
legend("Nonhierarchical","(\tau, \alpha)", "(\tau, M)", "(\alpha, M)", "(\tau, \alpha, M)");
xlabel('Sample number')
ylabel('Median accuracy')
hold off;
end