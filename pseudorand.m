k = ceil(log(n));
s = zeros(2,k);
t = linspace(0,n-1,n);
for i = 1:k
    c = randi(n - 1);
    cinv = modminv(c,n);
    p = mod(t*c,n)+1;
    permf = f(p);
    s(:,i)= mod(sparsefft(permf,2,n)*cinv,n);
end

s

