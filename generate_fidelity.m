function label_data = generate_fidelity(percent, truth, num_data)
%Fidelity (for MNIST data set)
%   Jumps specifies the changing points in the data

    y = randsample(num_data, floor(percent*num_data));
    label_data = zeros(num_data, 1);
    label_data(y) = truth(y);
end