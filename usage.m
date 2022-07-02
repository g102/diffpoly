% create an analytical signal as a random sum of three Fourier terms
x = linspace(0, 1, 1001);
y = zeros(size(x));
for jj = 1:3
	y = y + sin(2 * pi * randi(10) * x);
end

% add Gaussian noise
yPlusNoise = y + 0.05 * randn(size(x));

% compute derivatives
dy_clean = gradient(y, x(2)-x(1)); % estimator of correct derivative
dy_noisy = gradient(yPlusNoise, x(2)-x(1)); % naive approach with FD on noise
dy_diffpoly = diffpoly(x, yPlusNoise, 21, 2); % approach with diffpoly

%% plot for comparison
figure;
axes('NextPlot', 'add', 'box', 'on');
plot(x, dy_clean, 'displayname', 'Derivative of clean signal', 'linewidth', 1);
plot(x, dy_noisy, 'displayname', 'Estimation with finite differences', 'linewidth', 0.5);
plot(x, dy_diffpoly, 'displayname', 'Estimation with polynomial fitting', 'linewidth', 1);
xlabel('x');
ylabel('y(x)');
legend('location', 'no');
