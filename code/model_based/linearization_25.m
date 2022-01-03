function [matrix1, matrix2] = linearization_25(x, w, flag, dt)
%LINEARISATION_25 Used by kalmanFilter25sequential().
%   Linearizes ECG model at specified point of interest.
%   INPUT:  x   ->  System state, vector of length 25
%           w   ->  Process noise vector
%           flag->  if 0 linearisation of process function
%                   if 1 linearisation of observation function
%           dt  ->  discrete time step
%   OUTPUT: matrix1 ->  derivative of process or observation function with
%                       respect to state x
%           matrix2 ->  derivative of process or observation function with
%                       respect to process noise w
%
% Copyright 2019 Institute for Electrical Engineering in Medicine, 
% University of Luebeck
% Eike Petersen, Julia Sauer
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the 
% "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
% OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR 
% THE USE OR OTHER DEALINGS IN THE SOFTWARE.

phi = x(1);

alphai = x(5:11);
bi = x(12:18);
thetai = x(19:25);

omega = w(1);

if flag == 0    % linear approximation of process function    
    A = eye(25,25);
    B = eye(25,25);
    
    %dtetai = rem(phi - thetai, 2 * pi);
    dtetai = mod(phi - thetai + pi, 2 * pi) - pi;
    assert(all(-pi <= dtetai & dtetai <= pi));
    
    % d f / d phi
    A(2,1) = -dt * sum(alphai(1:2) * omega ./ (bi(1:2) .^2) .* (1 - dtetai(1:2) .^2 ./ (bi(1:2) .^2)) .* ...
        exp(-dtetai(1:2) .^2 ./ (2 * bi(1:2) .^2)));
    A(3,1) = -dt * sum(alphai(3:5) * omega ./ (bi(3:5) .^2) .* (1 - dtetai(3:5) .^2 ./ (bi(3:5) .^2)) .* ...
        exp(-dtetai(3:5) .^2 ./ (2 * bi(3:5) .^2)));
    A(4,1) = -dt * sum(alphai(6:7) * omega ./ (bi(6:7) .^2) .* (1 - dtetai(6:7) .^2 ./ (bi(6:7) .^2)) .* ...
        exp(-dtetai(6:7) .^2 ./ (2 * bi(6:7) .^2)));
    
    I = [2,27,53,78,103,129,154]; % repetitive matrix positions
    
    % d f / d alphai
    A(I + 100) = -dt * omega ./ (bi .^2) .* dtetai .* exp(-dtetai .^2 ./ (2 * bi .^2));
    
    % d f / d bi
    A(I + 275) = -2 * dt * alphai * omega ./ (bi .^3) .* exp(-dtetai .^2 ./ (2 * bi .^2)) .* ...
        (dtetai .^2 ./ (2 * bi .^2) - 1);
    
    % d f / d thetai
    A(I + 450) = -dt * alphai * omega ./ (bi .^2) .* exp(-dtetai .^2 ./ (2 * bi .^2)) .* (dtetai .^2 ./ (bi .^2) - 1);
    
    matrix1 = A;
    
    % d f / d omega
    B(1,1) = dt;
    B(2,1) = -dt * sum(alphai(1:2) ./ (bi(1:2) .^2) .* dtetai(1:2) .* exp(-dtetai(1:2) .^2 ./ (2 * bi(1:2) .^2)));
    B(3,1) = -dt * sum(alphai(3:5) ./ (bi(3:5) .^2) .* dtetai(3:5) .* exp(-dtetai(3:5) .^2 ./ (2 * bi(3:5) .^2)));
    B(4,1) = -dt * sum(alphai(6:7) ./ (bi(6:7) .^2) .* dtetai(6:7) .* exp(-dtetai(6:7) .^2 ./ (2 * bi(6:7) .^2)));
    
    matrix2 = B;
    
elseif flag == 1    %linear approximation of observation function
    
    M = eye(4,25);
    N = eye(4,4);
    
    matrix1 = M;
    matrix2 = N;
    
end
    

end

