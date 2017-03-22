function [ m ] = chirp( c,b,N )
%generates a chirp corresponding to the line L(a,b)
m= zeros(N,1);
c=mod(c,N);
b=mod(b,N);
if c~=0 
    for k=1:N
        m(k) = exp((2*pi*1i*((c*(k-1)^2*(N+1)/2)+b*(k-1)))/N)/(sqrt(N));
    end
elseif b ~=0
    for k=1:N
        m(k) = exp((2*pi*1i*(b*(k-1)))/N)/(sqrt(N));
    end
else
    m=ones(N,1)/(sqrt(N));
end
    
end

