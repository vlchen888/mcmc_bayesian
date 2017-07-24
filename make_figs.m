% To make figures after experiments
function make_figs(record, sigma_arr, fidelity_arr, color)
    for a = 1:length(fidelity_arr)
        subplot(1,3,a)
        % change me
        ylabel('average \tau')
        
        xlabel('sigma')
        title(sprintf('Percent fidelity = %.1f%%', fidelity_arr(a)*100))
        hold on
        curr_fidelity = fidelity_arr(a);
        to_plot = zeros(length(sigma_arr), 1);
        err_neg = zeros(length(sigma_arr), 1);
        err_pos = zeros(length(sigma_arr), 1);
        for b = 1:length(sigma_arr)
            curr_sigma = sigma_arr(b);

            TF1 = record(:,2)==curr_sigma;
            TF2 = record(:,3)==curr_fidelity;

            TFall = TF1 & TF2;

            trials = record(TFall, :);
            %p_mean = median(trials(:, 4));
            tau_mean = median(trials(:, 5));
            l = prctile(trials(:,5), 25);
            r = prctile(trials(:,5), 75);
            
            to_plot(b) = tau_mean;
            err_neg(b) = tau_mean - l;
            err_pos(b) = r - tau_mean;
        end
        errorbar(sigma_arr, to_plot, err_neg, err_pos, color, 'LineWidth', 2)
        axis([0 0.12 0 10]);
    end

end