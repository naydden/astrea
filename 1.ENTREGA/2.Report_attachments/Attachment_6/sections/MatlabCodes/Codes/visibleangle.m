function [Bo] = visibleangle (h,eo)

d=zeros(length(h),length(eo));
nsats_opt=d;
Re = 6371; 

for i = 1:length(h)
    for j = 1:length(eo)
        
        d(i,j)= Re * ( sqrt( ( (h(i)+Re)/Re )^2 - ( cosd(eo(j)) )^2 ) - sind(eo(j)) );
        Bo = asind( d(i,j)*cosd(eo(j)) / (h(i)+Re) );
        
    end
end

end
