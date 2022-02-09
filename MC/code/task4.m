% task 4

clear variables; close all;
javaaddpath('MCJava640.jar');

% mc=se.liu.imt.mcdata.VectorDataSimpleMode('task4_vessel_cyl.bin');
mc=se.liu.imt.mcdata.VectorDataSimpleMode('task4_tissue_cyl.bin');
X = mc.getXArray; % store the x-position in X
Z = mc.getZArray; % store the z-position in Z
Ne = mc.emitted; % store total n.o. emitted photons in Ne
Ntot = mc.detectedTotNo; % store total n.o. detected photons in Ntot

x_pos = [];
z_pos = [];

for i = 1:length(X)
    
    x = X(i);
    z = Z(i);

    x_pos = [x_pos x];
    z_pos = [z_pos z];
    
end

figure(1)
scatter(x_pos, z_pos, '.')
xlabel('x position')
ylabel('z position')

relativeNDetected = length(x_pos) / Ne;
