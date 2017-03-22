function [ b ] = modminv( a,n )
%compute the multiplicative inverse of a mod n, assuming a is invertible
%ie. gcd(a,n)=1
a = mod(a,n);
if a==0
    fprintf('this sucks')
    return
end
coeffs = eye(2); %coeffs in n and a (in that order)
r=100;
while r>0
    r = mod(n,a);
    q = floor(n/a);
    temp1= coeffs(1,1)-q*coeffs(2,1);
    temp2= coeffs(1,2)-q*coeffs(2,2);
    coeffs(1,:) = coeffs(2,:);
    coeffs(2,:) = [temp1, temp2];
    n=a;
    a=r;
end
    b=coeffs(1,2);
end

