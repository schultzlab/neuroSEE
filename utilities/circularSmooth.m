% Written by Simon Schultz
% smooth function that deals nicely with data being circular

function y = circularSmooth(x,w)

L = length(x);
postpend = x(1:w);
prepend = x(end-w+1:end);
x_widened = [prepend x postpend];
x_sm = smooth(x_widened,w);
y = x_sm(w+1:w+L);
