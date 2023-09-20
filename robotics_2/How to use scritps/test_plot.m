x = linspace(-10, 10, 100);
y1 = sin(x);
y2 = cos(x);

plot(x, y1)
hold on
%% 

plot(x, y2)
hold off

xlabel('x')
ylabel('y')
title('Plot of sin(x) and cos(x)')
legend('sin(x)', 'cos(x)')