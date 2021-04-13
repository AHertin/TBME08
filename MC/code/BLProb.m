function prob = BLProb(ma, ms, z)
% Probability function (Beer-Lambert)

    prob = -(1 / (ma + ms)) * log(1 - z);

end

