function [usvd1,usvd2,usvd3,usvd4] = usvdetection( USVP,Rusv_up )
%ODETECTION 此处显示有关此函数的摘要
%   此处显示详细说明

n = length(USVP(1,:));
usvd1 = [];usvd2 = [];usvd3 = [];usvd4 = [];
j=1;
for i=1:n
    if i~=1
        if norm(USVP(:,1)-USVP(:,i))<=Rusv_up
            usvd1(:,j) = USVP(:,i);
            j = j+1;
        end
    end
end
j=1;
for i=1:n
    if i~=2
        if norm(USVP(:,2)-USVP(:,i))<=Rusv_up
            usvd2(:,j) = USVP(:,i);
            j = j+1;
        end
    end
end
j=1;
for i=1:n
    if i~=3
        if norm(USVP(:,3)-USVP(:,i))<=Rusv_up
            usvd3(:,j) = USVP(:,i);
            j = j+1;
        end
    end
end
j=1;
for i=1:n
    if i~=4
        if norm(USVP(:,4)-USVP(:,i))<=Rusv_up
            usvd4(:,j) = USVP(:,i);
            j = j+1;
        end
    end
end

end

