function scatter_twomoons_classify(data, final_avg, label_data)
    colormap(redbluecmap(5))
    colors = -sign(final_avg); % to make blue = +, red = -
    scatter(data(:,1), data(:,2), 5 , colors)
    hold on
    for i = 1:length(label_data)
        if label_data(i) == 1
            scatter(data(i,1),data(i,2), 20, 'b', 'd', 'filled');
        elseif label_data(i) == -1
            scatter(data(i,1),data(i,2), 20, 'r', 'd', 'filled');
        end
    end
    hold off
end
