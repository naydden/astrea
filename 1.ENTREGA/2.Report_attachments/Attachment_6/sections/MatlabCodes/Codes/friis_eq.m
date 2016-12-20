function [ Pr,radius, FSL] = friis_eq( Gt, Gr, f, Pw, r_start, r_end, n,Lprop)
%FRIIS_EQ: Calculate the received on ground power in dBm.
%   Gt, Gr in dBi
%   f in Mhz
%   Pw in watts

%initiate
i=1;
Pr=zeros(size(r_start:n:r_end));
radius=zeros(size(Pr));
FSL=zeros(size(Pr));
P = 10*log10(Pw/1e-3);
for r=r_start:n:r_end
    
    Ls = 32.45+20*log10(r)+20*log10(f); %r-->Km; f-->Mhz%Spatial losses
    Pr(i) = P + Gt + Gr - Ls - Lprop;
    radius(i) = r;
    FSL(i)=Ls;
    %fprintf('Received power %fdBm \n',Pr(i));
    i=i+1;
end


end

