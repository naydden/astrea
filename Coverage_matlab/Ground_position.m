function [X] = Ground_position(lambda,long,t)

long=long*pi/180;
lambda=lambda*pi/180;
Rg=cos(lambda);
hg=sin(lambda);
w=2*pi/(24*60);
l=length(t);
hg=hg*ones(1,l);

X=[Rg*cos(w*t+long);Rg*sin(w*t+long);hg];

end

