function AF = AmbigFunc( f1 , f2 , c)
% Computes the ambiguity function for functions f1 and f2
% on Z/pZ
p = length(f1);
AF = zeros(p);
if c~=floor(c)
    sprintf('bad c value');
    return
end
if c>=0
    c = mod(c,p);
    for i=1:p
        Y=Ambigline(f1,f2,c,(i-1),p);
        for j=1:p
            if mod((j-1)*c+(i-1),p)== 0 
                AF(j,1) = Y(j); 
            else
                AF(j,mod((j-1)*c+(i-1),p)+1) = Y(j); 
            end
        end
    end
else
    for i=1:p
        Y=Ambigline(f1,f2,-1,(i-1),p)*1/sqrt(p);
        AF(i,:) = Y(:);
    end
end
AF=circshift(AF,[floor(p/2)-1,floor(p/2)-1]);
end


