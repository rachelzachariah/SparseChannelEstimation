function [dftk] = NonRecursiveSfft(N,ktemp,ak1,ak2,int)

A = zeros(log2(N),ceil(log2(log2(N))),3);
B = zeros(log2(N),ceil(log2(log2(N))),2);
Nvec = repmat(linspace(log2(N),1,log2(N))',[1,ceil(log2(log2(N)))]);
omega = exp(2*pi*1i./(2.^Nvec*2));
dftk = zeros(2*ktemp,2);

for i=1:2*ktemp   
    A(:,:,1) = ak1(:,:,i);
    A(:,:,2) = ak2(:,:,i);
    A(:,:,3) = int;    
    B(:,:,1) = A(:,:,1).*omega.^(int-1);
    B(:,:,2) = A(:,:,2).*omega.^(int-1+(2.^Nvec)/2);
    bits = zeros(1,log2(N));
    for k= 1:log2(N)
        m = log2(log2(N));
        if k==1 || (k>1 && bits(log2(N)+1-k+1)==0) 
            a = abs(A(k,:,1)-A(k,:,2)) - abs(A(k,:,1)+A(k,:,2));
        else
            a = abs(B(k,:,1)-B(k,:,2)) - abs(B(k,:,1)+B(k,:,2));
        end
        b = a<=0;
        if sum(b) <= floor(m/2)
            bits(log2(N)+1-k) = 0;
        else
            bits(log2(N)+1-k) = 1;
        end
    end
    dft2 = bin2dec(int2str(bits));
    coeff=0;
    omega = exp((2*pi*1i)/N);
    for j=1:log2(N)
        tempint = int(j,1)-1;
        coeff = coeff + ak1(j,i)*omega^(-tempint) +ak2(j,i)*omega^(-tempint-(N/2^N));
    end
    coeff = coeff/(log2(N)); %coefficient estimate
    dftk(i,:) = [coeff,mod(dft2,N)];
end   
end
    