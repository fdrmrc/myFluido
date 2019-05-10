function I=trapeziNonEq(xi,yi)

m=length(xi)-1;
I=0;
for i=1:m
    I=I+((xi(i+1)-xi(i))/2)*(yi(i+1)+yi(i));
 end
#Per calcolare questo si poteva anche usare il comando trapz presente in octave


