function L = compute_laplacian_standard(X, p, q, l)
    %Compute Laplacian matrix
    
    lap = 'un';
    
    W = exp(-pdist2(X, X, 'minkowski', p).^q/(2*l^2));
    D = diag(sum(W));
    
    L = D - W;
    if strcmp(lap, 'sym')
        L = sqrtm(inv(D)) * L * sqrtm(inv(D));
    elseif strcmp(lap, 'rw')
        L = inv(D) * L;
    end
end

function d = dist(x, y, p, q)
%Lp of x - y raised to the qth power
    d = norm(x - y, p)^q;
end

function w = weight(x, y, l, p, q)
%weight function with given metric and length-scale parameter
    w = exp(-dist(x, y, p, q)/ 2 / l^2);
end