function label_data = generate_voting_fidelity( number )
%Label some senators
%   percent = percent fidelity
y = randsample(435, number);
label_data = zeros(435, 1);
for i=1:length(y)
    if y(i) <= 267
        label_data(y(i)) = 1;
    else
        label_data(y(i)) = -1;
    end
end

end

