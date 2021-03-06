function mnist_heatmap(labels, truth, digs, txt)
    %Make MNIST heatmap
    nums_labels = labels * (1:length(digs))';
    nums_truth = truth * (1:length(digs))';
    
    TrueDigit = digs(nums_truth)';
    LabeledAs = digs(nums_labels)';
	count = zeros(length(TrueDigit),1);
    for i=1:length(digs)
        inds = find(TrueDigit==digs(i));
        count(inds) = 1/length(inds);
    end
    
    same = TrueDigit == LabeledAs;
    TrueDigit(same) = [];
    LabeledAs(same) = [];
    count(same) = [];
    

    
    
    t = table(LabeledAs,TrueDigit,count);
    hmo = heatmap(t, 'LabeledAs', 'TrueDigit','ColorVariable','count','ColorMethod','sum', 'ColorMap', parula);
    title(txt);
end

