function [ majorityvote ] = majorityVote( akj )
%does a majority vote on for the bit in question
%     counter = zeros(2);
%     counter(2,1) = 1;
    sz = size(akj);
    m = sz(1,2);
    
    a= abs(akj(1,:)-akj(2,:))-abs(akj(1,:)+akj(2,:));
    
    b = a<=0;
    if sum(b) <= floor(m/2)
        majorityvote = 0;
    else
        majorityvote = 1;
    end
    
%     for i = 1:m
%         if abs(akj(1,i)-akj(2,i)) <= abs(akj(1,i)+akj(2,i))
%             counter(2,2)=counter(2,2)+1;
%         else counter(1,2) = counter(1,2)+1;
%         end
%     end
%     if counter(2,2)< counter(1,2)
%         majorityvote = 0;
%     else majorityvote = 1;
%     end
end

