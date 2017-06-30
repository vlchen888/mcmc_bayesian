function plotBar(avg)

    posx = zeros(1,1);
    posy = zeros(1,1);
    ppos = 1;
    negx = zeros(1,1);
    negy = zeros(1,1);
    pneg = 1;
    for i=1:length(avg)
        if avg(i) > 0
            posx(ppos) = i;
            posy(ppos) = avg(i);
            ppos = ppos + 1;
        else
            negx(pneg) = i;
            negy(pneg) = avg(i);
            pneg = pneg + 1;
        end
    end
    hold on
    bar(posx, posy, 'b')
    bar(negx, negy, 'r')
    hold off
end