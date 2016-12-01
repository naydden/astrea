function [L] = Links(Xs,Xg,e_min)
e_min=e_min*pi/180;
l=size(Xg,2);
N_sat=size(Xs,2);
alpha_lim=pi/2-e_min;

L=zeros(1,l);

for i=1:l
    for j=1:N_sat
        Xgs=Xs(:,j,i)-Xg(:,i);
        alpha(j,i)=acos((Xgs'*Xg(:,i))/(norm(Xgs)));
        if alpha(j,i)<alpha_lim
            L(i)=L(i)+1;
        end
    end
end
end


