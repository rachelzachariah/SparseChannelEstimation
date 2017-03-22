function [ tfshift ] = TFshift( f, tau , omega )
%Takes input and returns time-frequency shift of w by
% (tau, omega)
N = length(f);
k = linspace(0,N-1,N)';
tfshift = exp((2*pi*1i*-omega*k)/N).*circshift(f,[-tau,0]);%e(-omega t)f(t - tau)
end

