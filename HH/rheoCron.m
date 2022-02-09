strdur = [0.05 131; 0.1 66; 0.2 33; 0.4 17; 0.7 10; 1.0 7];

s = strdur';

hold on
plot(s(1,:), s(2,:))
axis([0.05 1 0 135]);
xline(0.5285, '-', {'Chronaxy = 0.528 ms'})
yline(7, '-', {'Rheobase = 7 \muA'})
yline(14, '-', {'2x Rheobase = 14 \muA'})
