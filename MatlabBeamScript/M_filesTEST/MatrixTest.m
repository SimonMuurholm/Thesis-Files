clc
clear all

K = [5 1 9; 4 5 6; 7 8 10]
phi = [1;2;3]
F = [5 2 4 9 8 7 6; 5 2 4 9 8 7 6; 5 2 4 9 8 7 6;]

a = phi'*F
b = phi'*K*phi

a/b