function [ tfshift,m,l,l0 ] = EstimChan(H,L,M,M0,a1,a2,a3,r,sparse )
%Returns the channel parameters using incidence 
%and cross methods.
%H is time-frequency shifted echo recieved
%S1,S2,S3 represent slopes of the lines L,M,M0 as required.

[tfshift, m, l, l0] = tripleInt(a1,a2,a3,r,H,L,M,M0,sparse);
 
end

function [ tfshift, m, l, l0 ] = tripleInt(a1,a2,a3,r,H,RL,RM,RM0,sparse)
    N = length(H);
    tfshift= zeros(0,0);
    L= chirp(a1,0,N);
    M= chirp(a2,0,N);
    M0 = chirp(a3,0,N);
    if sparse==0
        tic
        [~,sortingl]= sort(Ambigline(L,RL,a2,0,N),'descend');
        m=sortingl(1:r)-1;
        [~,sortingm]= sort(Ambigline(M,RM,a1,0,N),'descend');
        l=sortingm(1:r)-1;
        [~,sortingl0]= sort(Ambigline(M0,RM0,a1,0,N),'descend');
        l0=mod((mod(sortingl0(1:r),N)-1)*(a1-a3),N);
        toc
    else
        tic
        m =  mod(sparseAmbigline(sqrt(N)*L,RL,r,N)*modminv(a2-a1,N),N);
        l =  mod(sparseAmbigline(sqrt(N)*M,RM,r,N)*modminv(a1-a2,N),N);
        l0 = mod(sparseAmbigline(sqrt(N)*M0,RM0,r,N),N);
        toc
    end

     for i=1:min(r,length(m))
         for j=1:min(r,length(l))
             c= mod((a2)*m(i)+a1*l(j)- a3*(m(i)+l(j)),N); % eqn of shifted line is y = y0+a3*(x-x0) so a point is on M0 
                                                % iff y coord = a2*x+
                                                % a3*(xcoord - x) 
                                                % for some x in l0.
              mod(l0-c,N);
%             mod([m(i)+l(j) (a2-N)*m(i)+a1*l(j)],N)
             if ~isempty(find(mod(l0-c,N)<500,1)) 
                x= mod(m(i)+l(j),N);
                y= mod(a2*m(i)+a1*l(j),N);
                tfshift = [tfshift',[x,y]']';
             end
         end 
     end
end

% function [ tfshift ] = doubleInt(a1,a2,r,H)
%     N = length(H);
%     tfshift= zeros(0,0);
%     L= chirp(a1,0,N);
%     M= chirp(a2,0,N);
%     [Am,sortingl]= sort(Ambigline(L,H,a2,0,N),'descend');
%     m=sortingl(1:r)-1
%     [Al,sortingm]= sort(Ambigline(M,H,a1,0,N),'descend');
%     l=sortingm(1:r)-1
%     hypothesis = @(j,k) Am(j)*exp(2*pi*1i*l(k)/N)- Al(k)*exp(2*pi*1i*m(j)/N) ;
%     for i=1:r
%         for j=1:r
%             %if hypothesis(i,j)<1/sqrt(N)
%                 x= mod(l(i)+m(j),N)
%                 y= mod(a2*m(j)+a1*l(i),N)
%                 tfshift = [tfshift',[x,y]']';
%             %end
%         end
%     end
% end

