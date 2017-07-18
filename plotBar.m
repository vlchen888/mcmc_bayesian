function plotBar(plot_me)
    neg_x = find(plot_me < 0);
    
    pos_x = find(plot_me > 0);
    
    zero_x = find(plot_me == 0);
    if ~isempty(zero_x)
        bar(zero_x, zeros(1,length(zero_x)), 'k');
        hold on
    end
    
    if ~isempty(neg_x)
        bar(neg_x, plot_me(neg_x), 'r')
        hold on
    end
    
    if ~isempty(pos_x)
        bar(pos_x, plot_me(pos_x), 'b');
    end
    hold off

end