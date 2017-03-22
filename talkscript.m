n=2^17-1;
a1=mod(randi(n),n);
a2=mod(randi(n),n);
a3=mod(randi(n),n);
L = chirp(a1,0,n);
M = chirp(a2,0,n);
M0 = chirp(a3,0,n);
S1= L+M+M0;
shift=[ randi(n-1),randi(n-1)]
R1=7*TFshift(S1,shift(1,1),shift(1,2));
RL = 7*TFshift(L,shift(1,1),shift(1,2));
RM = 7*TFshift(M,shift(1,1),shift(1,2));
RM0 =7*TFshift(M0,shift(1,1),shift(1,2));
tic
tfshift = EstimChan(R1,RL, RM, RM0,a1,a2,a3,1,1,0)
toc
tic
tfshiftsparse = EstimChan(R1,RL, RM, RM0,a1,a2,a3,1,1,1)
toc