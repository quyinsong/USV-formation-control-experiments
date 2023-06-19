function [od1,od2,od3,od4] = odetection( USVP,OP,RO_up )
%ODETECTION 此处显示有关此函数的摘要
%   此处显示详细说明

n = length(OP(1,:));
od1 = [];od2 = [];od3 = [];od4 = [];
j=1;
for i=1:n
    if norm(USVP(:,1)-OP(:,i))<=RO_up(i)
        od1(j,:) = [OP(:,i)',0,0,0,0];
        j = j+1;
    end
end
j=1;
for i=1:n
    if norm(USVP(:,2)-OP(:,i))<=RO_up(i)
        od2(j,:) = [OP(:,i)',0,0,0,0];
        j = j+1;
    end
end
j=1;
for i=1:n
    if norm(USVP(:,3)-OP(:,i))<=RO_up(i)
        od3(j,:) = [OP(:,i)',0,0,0,0];
        j = j+1;
    end
end
j=1;
for i=1:n
    if norm(USVP(:,4)-OP(:,i))<=RO_up(i)
        od4(j,:) = [OP(:,i)',0,0,0,0];
        j = j+1;
    end
end

end

