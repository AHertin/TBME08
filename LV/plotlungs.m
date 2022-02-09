close all;

plot(konc(rl,:))
hold on,
plot(konc(ll,:))
legend("right", "left", "location", "best")