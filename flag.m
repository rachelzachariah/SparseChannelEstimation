function [t,w] = flag(H,S1,S2,S3,r )
%Computes the flag lines and return slope and 
%intercept
t=0.2; %set threshold
p = length(H);
if S2.equals('none')
   l1 = Ambigline(S1,H,-1,0); 
   [~,l1]=max(l1, t);
   l1= sort(max(l1, t),'descend');
   l1 = l1(1:r);
   l2 = Ambigline(S1,H,-1,p-1); 
   l2= sort(max(l2, t),'descend');
   l2 = l2(1:r);
   if sum(any(l1))<sum(any(l2))
       w=l2;
   else w=l1;
   end
    
elseif S3.equals('none')
    
else
end

end

