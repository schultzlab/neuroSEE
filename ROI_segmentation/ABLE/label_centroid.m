function label_centroid(stat,col)
    count = 0;
    for x = 1: numel(stat)
    plot(stat(x).Centroid(1),stat(x).Centroid(2),col,'LineWidth',1);
    
    h = text(stat(x).Centroid(1),stat(x).Centroid(2),num2str(x));
    set(h,'Color',col(1));
    count = x;
    end


    
end

