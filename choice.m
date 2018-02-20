function y=choice(x)
for i=1:length(x)
    if x(i)>0 && x(i)<1
        y=x(i);
    end
    
   if exist('y')==0
       y=0;
   end
end
end