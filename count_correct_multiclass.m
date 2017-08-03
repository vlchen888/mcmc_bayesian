function p = count_correct_multiclass(final_labels, label_data, truth)

remainder = length(truth) - sum(sum(abs(label_data)));

wrong = sum(abs(final_labels-truth), 2)/2;
wrong = wrong(sum(label_data, 2)==0);
wrong = sum(wrong);
p = (remainder - wrong) / remainder;

end
