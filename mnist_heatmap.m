function mnist_heatmap( labels, truth, digs)
    %Make MNIST heatmap

    [nums_labels, ~] = find(labels);
    [nums_truth, ~] = find(truth);
    
    Digit = digs(nums_truth)';
    LabeledAs = digs(nums_labels)';
    
    
    t = table(LabeledAs, Digit);
    heatmap(t, 'LabeledAs', 'Digit');

end

