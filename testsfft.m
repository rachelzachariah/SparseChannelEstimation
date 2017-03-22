n = 2^17;
ktemp=3;
k=1;

%create signal:

% option 1
a=zeros(n,1);
freq1 = randi(n)-1;
freq2 = randi(n)-1;

for j=1:n
    a(j) = exp(2*pi*1i*(j-1)*(freq1)/n)+exp(2*pi*1i*(j-1)*(freq2)/n);%+exp(2*pi*1i*(j-1)*160/n);
end

temp = sparsefft(a,2,n)

tic
[~,in] = sort(abs(fft(a)));
in(1:2);
toc
%mod(temp(:,2),n)
    mod(freq1,n)
    mod(freq2,n)

%option2
%a=Ltemp.*(conj(Rtemp));

%create filter:
% 
% filter=zeros(n,1);
% for i=1:n %generate the significant coeffs of filter roughly |supp|=k
% if i < floor(k/2.0)+1
% filter(i)= exp(-(n/k)*(pi*(i-1)^2)/n);
% elseif i> n-floor(k/2.0)
% filter(i)= exp(-(n/k)*(pi*(n-i)^2)/n);
% end
% end

%create copies fo shifted filter:

% shiftfilter=zeros(n,k);
% shiftfilter(:,1)=filter;
% filterf= zeros(n,k);
% for i=1:k
% %     filterf(:,i)=cconv(a,shiftfilter,n);
% %     sparsefft(filterf,1,n);
%     for j=1:floor(k/2.0)+1
%         shiftfilter(j,i+1) = shiftfilter(j,i)*exp(2*pi*1i*(j-1)*(1/k));
%     end
%     for j=n-floor(k/2.0):n
%         shiftfilter(j,i+1) = shiftfilter(j,i)*exp(2*pi*1i*(j-1)*(1/k));
%     end 
% end

% ak1 = zeros(n,k);
% ak2 = zeros(n,k);
% atemp=zeros(n,1);
% for j=0:n-1 %generate samples of filtered versions of a each with one alive frequency for bit by bit
%         int1 = j; 
%         h2 = nonzeros(circshift(a,[-j,0]).*filter);
%         h1=zeros(k,1);
%         for i=0:k-1 
%               if i+1 <= floor(k/2.0)
%                 shiftint1 = mod(i+int1,n);
%                 h1(i+1) = a(shiftint1+1)*filter(i+1);
%               else
%                 shiftint1 = mod(n-(k-i)+int1,n);
%                 h1(i+1) = a(shiftint1+1)*filter(n-(k-i)+1);
%               end
%         end
%         ak1(j+1,:)= fft(h1,k); %creates samples for the n-th bit for each of the k filtered versions
%         ak2(j+1,:) = fft(h2,k);        
% end
% for i=1:k
%     sparsefft(ak1(:,i),1,n)
%     sparsefft(ak2(:,i),1,n)
% end

% ak1temp=zeros(log2(n),3);
% for i=1:k
% for j=1:log2(n)
%     int1 = randi(n);
%     int2 = mod(int1-1+(n/2^j),n)+1;
%     ak1temp(j,:) = [ak1(int1,i),ak1(int2,i),int1];
% end
% mod(sparsefftinternal(ak1temp,1,0,n),n)
% end

% ak1 = zeros(log2(n),k);
% ak2 = zeros(log2(n),k);
% int = zeros(log2(n),k);
% for j=1:log2(n) %generate samples of filtered versions of a each with one alive frequency for bit by bit
%         int1 = randi(n)-1; 
%         int2 = mod(int1+(n/2^j),n);
%         h1=zeros(k,1);
%         h2 = zeros(k,1);
%         for i=0:k-1 
%               if i+1 <= floor(k/2.0)
%                 shiftint1 = mod(i+int1,n);
%                 shiftint2 = mod(i+int2,n);
%                 h1(i+1) = a(shiftint1+1)*filter(i+1);
%                 h2(i+1) = a(shiftint2+1)*filter(i+1);
%               else
%                 shiftint1 = mod(n-(k-i)+int1,n);
%                 shiftint2 = mod(n-(k-i)+int2,n);
%                 h1(i+1) = a(shiftint1+1)*filter(n-(k-i)+1);
%                 h2(i+1) = a(shiftint2+1)*filter(n-(k-i)+1);
%               end
%         end
%         int(j)=int1+1;
%         ak1(j,:)= fft(h1,k); %creates samples for the n-th bit for each of the k filtered versions
%         ak2(j,:) = fft(h2,k);        
% end
% 
% for i=1:k
% mod(sparsefftinternal([ak1(:,i) ak2(:,i) int],0,n),n)
% end


% attempt at majority vote counter:

% rtemp = 1;%2^ceil(log2(k));
% m = zeros(rtemp*log2(n),2); %counter
%     for i=1:log2(n) 
%         mtemp=sparsefftinternal([ak1 ak2 int],rtemp,0,n);
%         [row,~]= find(m(:,2));  
%         for j =1:rtemp
%             check = sum(m(:,2)==mtemp(2,j));
%             max = max(row);
%             if check == 0                   %new frequency 
%                 m(max+1,:) = [mtemp(2,j),1];%initialize count of new freq  
%             else  
%                 for k=1:max
%                     if m(k,1)== mtemp(2,j)
%                         m(k,2)=m(k,2)+1;  %increase count of existing freq
%                     end
%                 end
%             end
%         end
%     end
