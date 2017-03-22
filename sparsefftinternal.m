function [dftk]=sparsefftinternal(ak,shift,N)
% ak a vector of length N=2^n, frequency-shifted by "shift", k the sparsity in Fourier domain
omega = exp(2*pi*1i*shift/(N*2));
int1 = ak(1,:,3)-1;
int2 = int1+N/2;
a1 = ak(1,:,1).*(omega.^int1);
a2 = ak(1,:,2).*(omega.^int2);
    if majorityVote([a1;a2]) %abs(a1-a2)<= abs(a1+a2)
         %'evensh'
         if N==2 
             dftk=2;
         else dftk= 2*sparsefftinternal(ak(2:log2(N),:,:),0,N/2);
         end  
     else 
         %'oddsh'
         if N==2
             dftk=1;
         else dftk= 2*sparsefftinternal(ak(2:log2(N),:,:),1,N/2)-1;
         end
    end
end

