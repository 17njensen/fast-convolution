rect=@(x,a) ones(1,numel(x)).*(abs(x)<a/2) % a is the width of the pulse
x=-0:0.001:1;
z = 0.1088;
y=rect(x,z*2);   
plot(x,y)
ylim([0 1.2]);
xlabel('Frequency(xPI rad/sample)', 'FontSize', 10);
ylabel('Amplitude','FontSize', 10);
title('Ideal Filter','FontSize', 10);