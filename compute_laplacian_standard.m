function [L, L_sym, L_rw] = compute_laplacian_standard(X, p, q, l)
    %compute Laplacian matrices
    [n, ~] = size(X);
    
    W = zeros(n);
    D = zeros(n);
    for i=1:n
        sum = 0;
        for j=1:n
            x = X(i,:);
            y = X(j,:);
            W(i, j) = weight(x, y, l, p, q);
            sum = sum + W(i, j);
        end
        D(i, i) = sum;
    end
    L = D - W; 
    L_sym = sqrtm(inv(D)) * L * sqrtm(inv(D));
    L_rw = inv(D) * L;
   %plotEigenvectors(L)
end

function d = dist(x, y, p, q)
%Lp of x - y raised to the qth power
    d = norm(x - y, p)^q;
end

function w = weight(x, y, l, p, q)
%weight function with given metric and length-scale parameter
    w = exp(-dist(x, y, p, q)/ 2 / l^2);
end