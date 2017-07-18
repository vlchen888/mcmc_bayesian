function plot_senator_traces(u_all, plot_me)
%PLOT SOME SENATORS
    
    width = length(plot_me)/2;
    for i = 1:length(plot_me)
        subplot(2, width,i)
        plot(movmean(sign(u_all(plot_me(i),:)), [length(u_all(plot_me(i),:)), 0]))
        ylabel(sprintf('Senator %d',plot_me(i)))
    end
end

