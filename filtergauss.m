N=97.0;
x=zeros(N,1);
r=sqrt(N);
num = linspace(-ceil(N/2),floor(N/2),N);
for i=1:N
    if mod(i-1,N) < N/2
        x(i,1) = exp(-(N/(r^2)*(i-1)^2)/N);
    else
        x(i,1) = exp(-(N/(r^2)*(N-i+1)^2)/N);
    end
end

% for i =1:N
%     x(i,1) = exp(-r*(pi*num(i)^2));
% end


% x=chirp(1,1,n);
x = circshift(x*(1/norm(x)),[ceil(N/2),0]);
 Z=abs(AmbigFunc(x,x,-1));
plot(num, abs(x));
hold on
y = abs(fft(x));
y = circshift(y * (1/norm(y)),[ceil(N/2),0]);
plot(num,y);
[X,Y] = meshgrid(-n/2+1:n/2,-n/2+1:n/2);
surf(X,Y,Z);
%  
