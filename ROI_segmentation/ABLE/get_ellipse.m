function [x1 y1]=get_ellipse(x0,y0,a,b,k)
%Draw the ellipse on a figure
%x0,y0 are the center of the ellipse;
%a is the major axis;b is the minor axis;
%k is the angle between the axis x and the major axis.
%Rotate transform matrix is [cos(k)  sin(k);
%                      -sin(k)  cis(k)]
sita=0:pi/20:2*pi;
x1=a*cos(sita)*cos(k)+b*sin(sita)*sin(k)+x0;
y1=-a*cos(sita)*sin(k)+b*sin(sita)*cos(k)+y0;
%figure,plot(x1,y1);
