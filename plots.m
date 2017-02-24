k = 
figure
subplot(4,1,1)
plot(peaks.peak13)
subplot(4,1,2)
plot(locs.loc13)
subplot(4,1,3)
plot(w.w13)
subplot(4,1,4)
plot(p.p13)


figure
subplot(4,1,1)
plot(peaks2.peak13)
subplot(4,1,2)
plot(locs2.loc13)
subplot(4,1,3)
plot(w2.w13)
subplot(4,1,4)
plot(p2.p13)

figure
for k = 1:6
subplot(6,1,k)
plot(fnorm.(sprintf('f%d',k)))
end

%randomqoms
figure
for k = 1:4
subplot(4,2,k)
list = 1:756;
num = [datasample(list,4)];
plot(allwins.(sprintf('win%d',num(k))))
title(['Quantity of Motion plot for',num2str(num(k))])
end


%random3dplots
figure
for k = 1:4
    subplot(4,2,k)
    list = 1:756;
    num = [datasample(list,4)];
    plot3(allwins.(sprintf('win%d',num(k)))(:,1),allwins.(sprintf('win%d',num(k)))(:,2),allwins.(sprintf('win%d',num(k)))(:,3))
    title(['motiontrace',num2str(num(k))])
end

%