% It seems in the implementation that certain choices of the ai are better 
% than others.Why??? EG n = 59, a1 = 53, a2=49, a3 = 42

% a1=mod(randi(n),n);
% a2=mod(randi(n),n);
% a3=mod(randi(n),n);
% a1=0;
% a2=20;
% a3=1;
L = chirp(a1,0,n);
M = chirp(a2,0,n);
M0 = chirp(a3,0,n);
S1= L+M+M0;
% S2= L+M;
% shift=round([ n/2 n/3 ; 0 n/4]);
% shift=[ randi(n-1),randi(n-1);randi(n-1),randi(n-1)];
shift
R1=7*TFshift(S1,shift(1,1),shift(1,2))+7*TFshift(S1,shift(2,1),shift(2,2));
RL = 7*TFshift(L,shift(1,1),shift(1,2))+7*TFshift(L,shift(2,1),shift(2,2));
RM = 7*TFshift(M,shift(1,1),shift(1,2))+7*TFshift(M,shift(2,1),shift(2,2));
RM0 =7*TFshift(M0,shift(1,1),shift(1,2))+7*TFshift(M0,shift(2,1),shift(2,2));
% R2=3*TFshift(S2,shift(1,1),shift(1,2))+5*TFshift(S1,shift(2,1),shift(2,2));

[tfshift,m1,l1,l10] = EstimChan(R1,RL, RM, RM0,a1,a2,a3,2,0);
%nonsptime(i) = 

[tfshiftsparse,m,l,l0] = EstimChan(R1,RL, RM, RM0,a1,a2,a3,2,1);
%sparsetime(i) = 

%tfshift2 = EstimChan(R2,a1,a2,a3,2,2)
 m2 = [mod((-a1*shift(1,1)+shift(1,2))*modminv(a2-a1,n),n), mod((a1*shift(2,1)-shift(2,2))*modminv(a1-a2,n),n)];
 l2 = [mod((-a2*shift(1,1)+shift(1,2))*modminv(a1-a2,n),n), mod((a2*shift(2,1)-shift(2,2))*modminv(a2-a1,n),n)];
 l02 = [mod(-a3*shift(1,1)+shift(1,2),n), mod(-a3*shift(2,1)+shift(2,2),n)];

tfshift
tfshiftsparse
% [m2 m']
% [l2 l']
% [l02 l0']
% hold on
% plot(abs(fft(L.*conj(RL))))
% plot(abs(fft(M.*conj(RM))))
% plot(abs(fft(M0.*conj(RM0))))
% hold off
% good vals for a2 :25977