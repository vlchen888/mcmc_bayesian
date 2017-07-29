function p = count_correct_multiclass(final_labels, label_data, truth)

remainder = length(truth) - sum(sum(abs(label_data)));

wrong = sum(abs(final_labels-truth))/2;
wrong = wrong(find(sum(label_data)==0));
wrong = sum(wrong);
p = (remainder - wrong) / remainder;

end
