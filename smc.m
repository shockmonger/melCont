load('mocoresults.mat')

%3d plots of handshapes for each window

    
%write adams typology for peaks, troughs of melody

%write parsons typology for peaks troughs of melody

%write fft feature vector output for peaks troughts of melody


%find similar melody - motion correspondances

[cidx2,cmeans2] = kmeans(allwins,2,'dist','sqeuclidean');


for i = 1:752
    LHZ(i,:) = allwins.(sprintf('win%d',i))(:,6)
    RHZ(i,:) = allwins.(sprintf('win%d',i))(:,3)
end


IDLHz = kmeans(LHZ, 16, 'distance', 'sqEuclidean')

figure
for i = 2:753
plot(normqom.(sprintf('n%d',i)),'g')
hold on
end
hold on
for i = 1:752
plot(LHZ(i,:),'b')
hold on
end
hold on
for i = 1:752
plot(RHZ(i,:),'r')
hold on
end

figure
subplot(2,1,1)
for i = 1:752
plot(LHZ(i,:))
hold on
end
subplot(2,1,2)
for i = 1:752
plot(RHZ(i,:))
hold on
end

figure 
for i = 1:16
    plot(p1.(sprintf('p%d',i)))
    hold on
end

p1 = {};
for i = 1:16
    p1.(sprintf('p%d',i)) = mat2gray(pitches.(sprintf('p%d',i)))
end

f2 = {};
for i = 1:16
   f2.(sprintf('p%d',i)) = interp(freqs.(sprintf('freq%d',i)),round(length(freqs.(sprintf('freq%d',i)))/20))
end

count = 1
for i = 2:753
    qomMat(:,count) = normqom.(sprintf('n%d',i))
    count = count+1;
end




[C, ia, ic] = unique(qomMat', 'rows')


for i = 2:753
    [U, I] = unique(qomMat(:,i), 'first');
end