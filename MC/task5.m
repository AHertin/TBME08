% task5

clear variables; close all;
javaaddpath('MCJava640.jar');

file1 = 'task5_1.bin';
file2 = 'task5_2.bin';
 
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
legend('with blood', 'only slab')

ratio = Id1./Id2;
x_ax = linspace(0,5,25);
figure(2)
plot(x_ax,ratio)
hold on;
yline(0.9)
title('Ratio between layer with blood and without')
ylabel('Ratio (%)')
xlabel('Distance from source (mm)')
