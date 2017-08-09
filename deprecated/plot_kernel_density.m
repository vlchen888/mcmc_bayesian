function plot_kernel_density()

    [tau_all, alpha_all, std] = modelB_learn_t_a();
    comp = ([tau_all ; alpha_all])';
    [bandwidth, density, X, Y] = kde2d(comp, 256, [0 0], [60 2.5]);
    figure(1)
    clf
    surf(X,Y,density,'LineStyle','none'), view([0,90])
    colormap hot, hold on, alpha(.3)
    set(gca, 'color', 'blue');
    xlabel('\tau');
    ylabel('\alpha');
    plot(comp(:,1),comp(:,2),'r.','MarkerSize',5)
    shg
end

