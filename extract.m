%Extract the values of geopotential height at grid points from the csv
%file (the fifth column), V is a row matrix with 10512 values
filename = '102412.csv';
M = csvread(filename,1,4)
V = M.'

%Create row matrix D s.t. the order is the same as the one used in 
%eniac.m
D = zeros(1,10512)

%Reorder the values (starts from the row y = 73, then y = 72, ..., 1)
for n = 1:73
    D(144*(n-1)+1:144*(n-1)+144) = V((73-n)*144+1:(73-n)*144+144)
end

%Save the matrix in asc format
dlmwrite('2012102412.asc',D)

