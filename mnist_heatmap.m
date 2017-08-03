function mnist_heatmap( labels, truth, digs)
    %Make MNIST heatmap
    nums_labels = labels * (1:length(digs))';
    nums_truth = truth * (1:length(digs))';
    
    TrueDigit = digs(nums_truth)';
    LabeledAs = digs(nums_labels)';
    count = zeros(length(truth),1);
    for i=1:length(digs)
        inds = find(TrueDigit==digs(i));
        count(inds) = 1/length(inds);
    end
    
    
    t = table(LabeledAs,TrueDigit,count);
    heatmap(t, 'LabeledAs', 'TrueDigit','ColorVariable','count','ColorMethod','sum'...
        , 'ColorMap', hot);

end

