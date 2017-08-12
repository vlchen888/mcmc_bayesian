function [ data, truth ] = generate_mnist_data(digs)
    K = length(digs);
    num = 4000;
    N0 = num*10/K;
    [imgs,labels] = readMNIST('train-images-idx3-ubyte','train-labels-idx1-ubyte',N0,0);

    data0 = reshape(imgs,400,N0)';
    data = [];
    for k=1:K
        data = [data;data0(labels==digs(k),:)];
    end

    N = size(data,1);

    %
    % Project onto first 50 pca components
    coef = pca(data);
    C = coef(:,1:50);
    data = (C*((C'*C)\(C'*data')))';
    %

    jumps = [0,cumsum(sum(labels==digs))];
    truth = zeros(N,K);
    for k=1:K
        truth(jumps(k)+1:jumps(k+1),k) = 1;
    end
end

