function p = test_fiedler_clustering_accuracy(sigma)
%Compute Fiedler Classification Accuracy
%   p = compute_fiedler_clustering_accuracy() computes the accuracy
%   Uses kmeans to identity 2 clusters from the Fiedler vector of the
%   Laplacian (currently unnormalized).

N = 2000;
data = moondata(1,100,N,sigma);

params = containers.Map;
params('p') = 2;
params('q') = 2;
params('l') = 1;
params('data') = data;
params('laplacian') = string('self tuning');


clustering = fiedler_clustering(params);
scatter_twomoons_classify(data, clustering, [])
actual = [zeros(floor(N/2)+1,1) - 1; zeros(N-(floor(N/2)+1),1) + 1];
p = count_correct(clustering, zeros(N,1), actual);
end

