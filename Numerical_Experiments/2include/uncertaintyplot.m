function h = uncertaintyplot(x,data,linestyle,color1,color2,ylog,s)

    y2 = [x, fliplr(x)];
    
    lowerbound = mean(data)-s*std(data);
    for i = 1:numel(lowerbound)
        if lowerbound(i) <= 0 && ylog == true
            lowerbound(i) = 0.01;
        end
    end
    
    inBetween = [(mean(data)+s*std(data))'; fliplr(lowerbound)'];
    fill(y2, inBetween,color2,'LineStyle','none'); hold on;
    h = plot(x,mean(data),'LineStyle',linestyle,'Color',color1,'LineWidth',1);
    if ylog == true
        set(gca, 'YScale', 'log')
    end
end