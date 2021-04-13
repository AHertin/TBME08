% task 4

clear variables; close all;
javaaddpath('MCJava640.jar');

mc=se.liu.imt.mcdata.VectorDataSimpleMode('task4.bin');
%mc=se.liu.imt.mcdata.VectorDataSimpleMode('task4_tissue_only.bin');
X = mc.getXArray; % store the x-position in X
Y = mc.getYArray; % store the y-position in Y
Z = mc.getZArray; % store the z-position in Z
Ne = mc.emitted; % store total n.o. emitted photons in Ne
Ntot = mc.detectedTotNo; % store total n.o. detected photons in Ntot
Na = mc.absorbed; % store total n.o. absorbed photons in Na

lg = length(X);
x_pos = [];
z_pos = [];

for i = 1:length(Y)
    
    x = X(i);
    y = Y(i);
    z = Z(i);
    
    if (0 < y && y < 100)
        if(0.95 <= z && z <= 1.05)
            if(sqrt(x^2 + (z-1)^2) < 0.05)
                x_pos = [x_pos x];
                z_pos = [z_pos z];
            end
        end
    end 
end

figure(1)
scatter(x_pos(1:1000), z_pos(1:1000), '.')
xlabel('x position')
ylabel('z position')

relativeNDetected = length(x_pos) / Ne;
