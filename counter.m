function [ cntr ] = counter(r,sparsity,error)
%majority vote of frequencies
r = abs(r);
r(:,2)=round(r(:,2));
[n,~]= size(r);
cntr = zeros(n,3);
numel=1;
for i=1:n
    found=0;
    for j=1:numel
        if cntr(j,3)== r(i,2);
            cntr(j,2)= cntr(j,2)+1;
            cntr(j,1) = cntr(j,1)+abs(r(i,1));
            found=1;
            break
        end
    end
    if found==0
        cntr(numel,:)=[r(i,1), 1, r(i,2)];
        numel=numel+1;
    end
end
cntr = cntr(1:numel-1,:);
cntr = sortrows(cntr,[-2,-1]);
% cntr(:,2)
% cntr(:,1)
% cntr(:,3)
[~,I]=sort(cntr(:,2),'descend');
cntr = cntr(I,:);
numfreqs=1; % keeps count of how many new frequencies have been detected
freqs = zeros(sparsity,3);
freqs(1,:) = cntr(1,:);
for i=1:numel-1
     found=0;
     for j = 1:sparsity %search for repeats
        if abs(freqs(j,3)-cntr(i,3))<= error
            freqs(j,2) = cntr(i,2)+1;
            freqs(j,1) = cntr(i,1)+freqs(j,1);
            found=1;
            break
        end
    end
    if found==0 %if not found start a new frequency counter
        numfreqs = numfreqs+1;
        freqs(numfreqs,:) = cntr(i,:);    
    end
    if numfreqs>=sparsity 
        break
    end  
end
cntr=freqs(1:min(sparsity,numfreqs),:,:);
end

