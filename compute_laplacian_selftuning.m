function L = compute_laplacian_selftuning(data)

lap = 'sym';                 % Laplacian: 'un','sym' or 'rw'
dist_type = 'euclidean';    % Distance for weights
q = 2;                      % Power of distance for weights
ell = 1;                    % Length-scale for weights
K = 10;                     % Number of nearest neighbours

N = size(data,1);
W=zeros(N,N); 

[IDX,DX] = knnsearch(data,data,'K',K,'Distance',dist_type,'IncludeTies',true);

for i=1:N
    for j=1:length(IDX{i})
        l = IDX{i}(j);
        W(i,l) = exp(-DX{i}(j)^q/(2*ell^q*DX{i}(K)*DX{j}(K)));
        %W(i,l) = exp(-DX{i}(j)^q/(2*ell^2));
    end
end
% Take care of 0/0 cases and symmetrize, just in case
W(isnan(W)) = 1;
W = (W+W')/2;

D=diag(sum(W));
DD=diag(1./sqrt(sum(W)));
L = D-W;

if strcmp(lap,'sym')
    L = DD*L*DD;
elseif strcmp(lap,'rw')
    L = DD*DD*L;
end
L = (L+L')/2;

end
