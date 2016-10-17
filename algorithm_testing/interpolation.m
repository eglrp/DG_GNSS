load('sat.mat');
load('signal.mat');
load('boundary.mat');

x(1:4) = 3600*date_ref(1:4,4) + 60*date_ref(1:4,5);
x = x';

for i = 1:10
    for j = 1:4
        y1(j, i) = XS_ref( i, 1, j);
        y2(j, i) = XS_ref( i, 2, j);
        y3(j, i) = XS_ref( i, 3, j);
    end
end

for i = 1:size(date,1)
    xi(i,1) = 3600*date(i,4) + 60*date(i,5) + date(i,6);
end

for i = 1:10
    yi_1(:,i) = pchip(x,y1(:,i),xi);
    yi_2(:,i) = pchip(x,y2(:,i),xi);
    yi_3(:,i) = pchip(x,y3(:,i),xi);
    dtS(i,:)  = pchip(x,dtS_ref(i,:),xi);
end

for i = 1:400
    XS(:,:,i) = [yi_1(i,:); yi_2(i,:); yi_3(i,:)];
end
