function label_data = generate_fidelity_multiclass(percent, truth, num_data, k)
%Fidelity multiclass (for MNIST data set)
%   Jumps specifies the changing points in the data

    y = randsample(num_data, floor(percent*num_data));
    label_data = zeros(num_data, k);
    label_data(y,:) = truth(y,:);
end