% task3

clear variables; close all;
javaaddpath('MCJava640.jar');

file1 = 'task3_5.bin';
file2 = 'task3_20.bin';

[Id1, x1] = calcIntensity(file1);
[Id2, x2] = calcIntensity(file2);

%%
figure(1)
semilogy(x1,Id1)
hold on;
semilogy(x2,Id2)
title('Scatter intensity by distance from source')
ylabel('Intensity (log)')
xlabel('Distance from source (mm)')
legend('\mu_s = 5', '\mu_s = 20')