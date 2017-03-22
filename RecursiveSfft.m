function [ dftk ] = RecursiveSfft( N,k,ak1,ak2,int)
% Internal function to perform the estimates the k most significant Fourier coefficients 
% given ak1, ak2

A = zeros(log2(N),ceil(log2(log2(N))),3);
dftk = zeros(2*k,2);

    for i=1:2*k
        A(:,:,1) = ak1(:,:,i);
        A(:,:,2) = ak2(:,:,i);
        A(:,:,3) = int;
        dft = sparsefftinternal(A,0,N);
        coeff=0;
        omega = exp((2*pi*1i)/N);
        for j=1:log2(N)
            tempint = int(j,1)-1;
            coeff = coeff + ak1(j,i)*omega^(-tempint) +ak2(j,i)*omega^(-tempint-(N/2^N));
        end
        coeff = coeff/(log2(N)); %coefficient estimate
        dftk(i,:) = [coeff,mod(dft,N)];
    end 

end

