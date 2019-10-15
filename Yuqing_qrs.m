% 
% Program info
% Program author: Yuqing Li, August 2019
% Program details: This matlab function generates a quadratic residue sequence.
% Inputs: 
% q: a prime number
% p: number of repitition
% Quadratic residue sequence: sn = (n^2)mod(q), where q is a prime number
% and (n^2)mod(q) is the least nonnegative remainder of modulo q. (Boundary and Medium Modelling Using Compact.....)
% The output y is a row vector.

function h = qrs(q,p)
if q<3 || ~isprime(q)
    error('qrd:primality', 'p must be an odd prime.');
else
    qrd = mod([0:q-1]'.^2,q);
end

h = qrd;
% repeat the QRD for p times
if p>1
    for k = 1:p
        h((k-1)*q+1:k*q) = h(1:q);
    end
end
end

% Reference:
% This function is modified from the work of Robert Dickson, available at https://uk.mathworks.com/matlabcentral/fileexchange/46890-qrseq-m