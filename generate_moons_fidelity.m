function label_data = generate_moons_fidelity(percent, N)
%Label moon data
%   percent = percent fidelity
%   N = number of data points

y = randsample(N, floor(percent*N));
label_data = zeros(N, 1);
for i=1:length(y)
    if (y(i)-1)*(2*pi/N) > pi
        label_data(y(i)) = 1;
    else
        label_data(y(i)) = -1;
    end
end

end

