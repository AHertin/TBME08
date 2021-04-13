%% Step size simulation

steps = 10^6;
muA   = 15;
muS   = 150;

Z = zeros(1,steps);

for i = 1:steps
    
    rnd = rand(1);
    
    Z(1,i) = BLProb(muA, muS, rnd);
    
end

%% Figure

histogram(Z,100)
title('Histogram of randomized step sizes')
xlabel('Step size (mm)')
ylabel('Number of occurences')

%% Comparison to theoretical

avg = mean(Z);
mfp = 1 / (muA + muS);

