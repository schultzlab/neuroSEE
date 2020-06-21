% plot raster of events with amplitude
% assumes event list X{c} for c cells
% X{c}.ampl, X{c}.time as vectors

function plot_amplraster(X,col,sz)
if nargin<2
    col = 'b';
end
if nargin<3
    sz=5;
end

C = length(X);
for i=1:C
    k = length(X{i}.time); % number of events for this cell
    yval = ones(1,k)*i;
    ev_size = X{i}.ampl./max(X{i}.ampl)*sz;
    scatter(X{i}.time,yval,ev_size,col,'filled'); 
    hold on
end