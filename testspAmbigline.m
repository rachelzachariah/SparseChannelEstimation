n=97;
a1=1;
L = chirp(1,0,n);
M = chirp(-1,0,n);
M0 = chirp(0,0,n);
S1= L+M+M0;
shift=[floor(0),floor(n/4);floor(3*n/4),floor(n/2)]
%shift=[ randi(n-1),randi(n-1);randi(n-1),randi(n-1)];
mod(shift(1,2)-shift(1,1),n) , mod(shift(2,2)-shift(2,1),n)
R1=1*TFshift(L,shift(1,1),shift(1,2))+1*TFshift(L,shift(2,1),shift(2,2));
R2 =TFshift(S1,shift(1,1),shift(1,2));%+TFshift(S1,shift(2,1),shift(2,2));
% tic
% round(sparseAmbigline(L*100,R2,1,n))
% toc
% tic
% [~,sortingl]= sort(Ambigline(L*100,R2,-1,0,n),'descend');
% m=sortingl(1:1)-1 
% toc
round(sparseAmbigline(L,R1,2,n))
% errorrate = 0;
% for i=1:100
%     if abs(131072-round(sparseAmbigline(L*100,R2,1,n)))>10
%         errorrate= errorrate+1;
%     end
% end
% errorrate= errorrate/100
% 

