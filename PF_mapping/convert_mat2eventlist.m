% Converts a cells x time matrix containing sparse events into a structure
% which lists for each cell the event times and amplitudes

function a = convert_mat2list(A)

C = size(A,1);
L = size(A,2);

for c=1:C
   idx = find(A(c,:)>0);
   a{c}.ampl = A(c,idx);
   a{c}.time = idx;
end