function [nRow, nCol] = getnRownCol(N)

fun = @(N,a,b) N > a && N <= b;
if N < 6
    nRow = 1; nCol = N;
elseif N == 6
    nRow = 2; nCol = 3;
elseif fun(N,6,8)
    nRow = 2; nCol = 4;
elseif fun(N,8,10)
    nRow = 2; nCol = 5;
elseif fun(N,10,12)
    nRow = 3; nCol = 4;
elseif fun(N,12,15)
    nRow = 3; nCol = 5;
elseif N == 16
    nRow = 4; nCol = 4;
elseif fun(N,16,20)
    nRow = 4; nCol = 5;
elseif fun(N,20,25)
    nRow = 5; nCol = 5;
elseif fun(N,25,30)
    nRow = 5; nCol = 6;
elseif fun(N,30,35)
    nRow = 5; nCol = 7;
elseif fun(N,35,40)
    nRow = 5; nCol = 8;
elseif fun(N,40,45)
    nRow = 5; nCol = 9;
elseif fun(N,45,50)
    nRow = 5; nCol = 10;
elseif fun(N,50,60)
    nRow = 5; nCol = 6;
elseif fun(N,60,70)
    nRow = 5; nCol = 7;
elseif fun(N,70,80)
    nRow = 5; nCol = 8;
elseif fun(N,90,100)
    nRow = 5; nCol = 10;
else
    nRow = 5; nCol = 8;
end
    