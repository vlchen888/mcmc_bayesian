
function S = compute_S_multiclass(u, k)
    [~, I] = max(u,[],2);
    A = eye(k);
    S = A(I, :);
end