function p = count_correct(final_avg, label_data, correct_labels)

remainder = length(correct_labels) - sum(abs(label_data));
cluster = sign(final_avg);
p = sum(abs(cluster-correct_labels).*(1-abs(label_data)))/2;
p = (remainder - p) / remainder;

end
