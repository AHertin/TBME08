function [Id, xax] = calcIntensity(filename)

    mc=se.liu.imt.mcdata.VectorDataSimpleMode(filename);
    X = mc.getXArray; % store the x-position in X
    Y = mc.getYArray; % store the y-position in Y
    Ne = mc.emitted; % store total n.o. emitted photons in Ne
    Ntot = mc.detectedTotNo; % store total n.o. detected photons in Ntot

    PhotonDists = (X.^2 + Y.^2);
    totMM = 5;
    distMM = 0.2;
    rings = totMM/distMM;
    Nd = zeros(1, rings);
    Id = zeros(1, rings);

    for j = 1:rings

        ring = j * distMM;

        for i = 1:Ntot
            pos = sqrt(PhotonDists(i));
            if (pos < ring)
                Nd(1, j) = Nd(1, j) + 1;
            end
        end

        if j > 1
            Nd(1, j) = Nd(1, j) - sum(Nd(1, 1:j-1));
        end

        Ad = pi * (ring^2 - (ring-distMM)^2);

        Id(1,j) = Nd(1,j) / (Ad * Ne);

    end

    %%
    xax = linspace(0,totMM,rings);

end

