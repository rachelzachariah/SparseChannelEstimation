function [ dftk ] = sparsefft( a, k ,N)
% a a vector of length N=2^n, k the sparsity in Fourier domain
if k>1
    ktemp = k;
    k = max(ktemp,8)/2;
end
if log2(N)-floor(log2(N))~=0
    Ntemp = N;
    N = 2^(ceil(log2(Ntemp))+1);
    a = [a;a;a;a];
    a = a(1:N);
else 
    Ntemp = N;
end
if k==1
    dftk = zeros(k,2);
    %only retain the entries of a that we need
    ak = zeros(log2(N),3,ceil(log2(log2(N))));
    
    for m = 1:ceil(log2(log2(N)))
        for j=1:log2(N)
            int1 = randi(N-1)+1;
            int2 = mod(int1-1+(N/2^j),N)+1;
            ak(j,:,m) = [a(int1),a(int2),int1];
        end
    end
    
    if abs(ak(1,1)-ak(1,2))<= abs(ak(1,1)+ak(1,2))
         %'even'
         if N==2 
             dftktemp=2;
         else dftktemp= 2*sparsefftinternal(ak(2:log2(N),:),0,N/2);
         end  
    else 
         %'odd'
         if N==2
             dftktemp=1;
         else dftktemp= 2*sparsefftinternal(ak(2:log2(N),:),1,N/2)-1;
         end
    end
    
    dftk(1,:)= [a(1),mod(dftktemp,N)]; %sparsefftinternal(ak,0,N);
    
else
    
    %pseudo random permutation:
     b = randi(N)-1;
     c = randi(N-1);
     cinv = modminv(c,N);
     
    %only retain the entries of a that we need
    dftk = zeros(2*k,2);
    filter=zeros(2*k,1);
    
    for i=1:k %generate the significant coeffs of filter roughly |supp|=k
        filter(i) = exp(-(N/(k^2))*(pi*(i-1)^2)/N);
    end
    j= k+1; 
    for i=N-k+1:N
        filter(j) = exp(-(N/(k^2))*(pi*(N-i)^2)/N);
        j=j+1;
    end
    
    ak1 = zeros(log2(N),ceil(log2(log2(N))),2*k);
    ak2 = zeros(log2(N),ceil(log2(log2(N))),2*k);
    int = zeros(log2(N),ceil(log2(log2(N))),1);

    for j=1:log2(N) %generate samples of filtered versions of a each with one alive frequency for bit by bit
       for l = 1:ceil(log2(log2(N)))
            int1 = randi(N)-1; 
            int2 = mod(int1+(N/2^j),N);
            h1=zeros(2*k,1);
            h2 = zeros(2*k,1);
            for i=0:2*k-1 
                  if i+1 <= k
                    shiftint1 = mod((i+int1),N);%mod(c*(i+int1),N);
                    shiftint2 = mod((i+int2),N);%mod(c*(i+int2),N);
                    h1(i+1) = a(shiftint1+1)*filter(i+1);%a(shiftint1+1)*filter(i+1);
                    h2(i+1) = a(shiftint2+1)*filter(i+1);%a(shiftint2+1)*filter(i+1);
                  else
                    shiftint1 = mod((N-(2*k-i)+int1),N);%mod(c*(N-(k-i)+int1),N);
                    shiftint2 = mod((N-(2*k-i)+int2),N);%mod(c*(N-(k-i)+int2),N);
                    h1(i+1) = a(shiftint1+1)*filter(i+1); %a(shiftint1+1)*filter(i+1); 
                    h2(i+1) = a(shiftint2+1)*filter(i+1); %a(shiftint2+1)*filter(i+1); 
                  end
            end
            int(j,l)=int1+1;
            ak1(j,l,:)= fft(h1,2*k); %creates samples for the n-th bit for each of the k filtered versions
            ak2(j,l,:) = fft(h2,2*k);      
       end
    end
    
    for i=1:2*k
        dft = sparsefftinternal([ak1(:,i) ak2(:,i) int],0,N)*Ntemp/N;
        coeff=0;
        omega = exp((2*pi*1i)/N);
        for j=1:log2(N)
            tempint = int(j,1)-1;
            coeff = coeff + ak1(j,i)*omega^(-tempint) +ak2(j,i)*omega^(-tempint-(N/2^N));
        end
        coeff = coeff/(log2(N)); %coefficient estimate
        dftk(i,:) = [coeff,mod(round(dft),N)];
    end   
    mtemp = counter(dftk,ktemp,10);
    m = mtemp(:,[1,3]);
    sz = size(mtemp);
    if ktemp>sz(1)
        fprintf('r is toooo big!');
        dftk = round(m(:,2));
    else
%       round(counter(mtemp,rtemp,30)*N/Ntemp);
        dftk = round(m(1:ktemp,2));
    end
    
end
end

function [dftk]=sparsefftinternal(ak,shift,N)
% ak a vector of length N=2^n, frequency-shifted by "shift", k the sparsity in Fourier domain
omega = exp(2*pi*1i*shift/(N*2));
int1 = ak(1,3)-1;
int2 = int1+N/2;
a1 = ak(1,1)*(omega^int1);
a2 = ak(1,2)*(omega^int2);
    if abs(a1-a2)<= abs(a1+a2)
         %'evensh'
         if N==2 
             dftk=2;
         else dftk= 2*sparsefftinternal(ak(2:log2(N),:),0,N/2);
         end  
     else 
         %'oddsh'
         if N==2
             dftk=1;
         else dftk= 2*sparsefftinternal(ak(2:log2(N),:),1,N/2)-1;
         end
    end
    
end
