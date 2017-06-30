function p = compute_fiedler_clustering_accuracy()
%Compute Fiedler Classification Accuracy
%   p = compute_fiedler_clustering_accuracy() computes the accuracy
%   Uses kmeans to identity 2 clusters from the Fiedler vector of the
%   Laplacian (currently unnormalized).

load('intertwine_moon.mat');
data = d;
[num_data, ~] = size(d);

params = containers.Map;
params('p') = 2;
params('q') = 2;
params('l') = 1;
params('data') = d;
params('laplacian') = string('un');


clustering = fiedler_clustering(params);
scatter_twomoons_classify(data, clustering, [])
actual = [zeros(1, num_data/2)-1, zeros(1,num_data/2)+1]';
wrong = min(norm(clustering - actual)^2/4, norm(clustering + actual)^2/4);
p = (num_data-wrong)/num_data;
end

