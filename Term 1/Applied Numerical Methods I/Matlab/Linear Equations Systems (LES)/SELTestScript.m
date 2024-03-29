format shortG
clear,clc
% A = [4 -1 0 0; -1 4 -1 0; 0 -1 4 -1; 0 0 -1 4];
% b = [1 1 1 1];
% x0 = [0 0 0 0];
% A = [10 -1 2 0; -1 11 -1 3; 2 -1 10 -1; 0 3 -1 8];
% b = [6 25 -11 15];
% x0 = [0 0 0 0];
% A = [2 1 -1; -3 -1 2; -2 1 2];
% b = [8 -11 -3];
% x0 = [0 0 0];
A = [1 0 0 1; 0 2 0 1; 0 0 3 1; 1 1 1 4];
b = [2 3 4 7];
x0 = [0 0 0 0];
w = 1.3;

x_jacobi = Jacobi(A, b, x0, 1e-7, 50);
x_jacobi

x_GS = GaussSeidel(A, b, x0, 1e-7, 50);
x_GS

x_SOR = SOR(A, b, x0, w, 1e-7, 50);
x_SOR