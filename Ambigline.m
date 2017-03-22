function AFi = Ambigline(f1,f2, c, w ,N)
if c>=0
    m1 = (sqrt(N)*chirp(N-c,N-w,N)).*f1;
    
    m2 = circshift(flipud((sqrt(N)*chirp(c,0,N)).*conj(f2)),[1,0]);    
    %the circshift is so that in flipping to get f_ the zero-th entry does not move.
    
    AFi= cconv(m1,m2,N);
else
    m2 = conj(circshift(f2,[w,0]));
    AFi = sqrt(N)*fft(f1.*m2);
end
end


