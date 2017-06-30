function [clusters] = fiedler_clustering(params)
    data = params('data');
    p = params('p');
    q = params('q');
    l = params('l');
    
    if params('laplacian') == string('self tuning')
        L = compute_laplacian_selftuning(data);
    elseif params('laplacian') == string('un')
        L = compute_laplacian_standard(data, p, q, l);
    end
    
    [phi, ~] = eig(L);
    clusters = kmeans(phi(:,2), 2)*2-3;
end