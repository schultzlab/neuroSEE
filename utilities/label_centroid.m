function label_centroid(stat,col)
    count = 0;
    for x = 1: numel(stat)
    plot(stat(x).Centroid(1),stat(x).Centroid(2),col,'LineWidth',1);
    
    h = text(stat(x).Centroid(1),stat(x).Centroid(2),num2str(x));
    set(h,'Color',col(1));
    count = x;
    end
    fprintf('\nTotal #of centroids %d\n',count);
end


% x_cen = zeros(200);
% y_cen = zeros(200);
% 
% for x = 1: numel(stat_r)
%    x_cen(x) = stat_r(x).Centroid(1);
%    y_cen(x) = stat_r(x).Centroid(2);
%     
% end
% 
% centroids = [x_cen,y_cen];