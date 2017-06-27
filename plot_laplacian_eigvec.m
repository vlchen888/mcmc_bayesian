
function plot_laplacian_eigvec(L)
    [V, ~] = eig(L);
    
    % Plot eigenvectors and determine clustering
    num_plots = 9;
    figure(1)
    for i = 1:num_plots
        subplot(3, 3, i);
        plot(1:length(V(:, i)), V(:, i))
    end
end
