x = zeros(100, 1);
ind=[12, 28, 44, 56, 67, 79, 33, 51, 90, 2];
for i=1:10
    x(ind(i)) = rand()*10;
end
h = [1, 2, 3, 4, 3, 2, 1]/16.' ;

y = conv(h, x) + 0.05*norm(x)*randn(106, 1);
A = convmtx(h, 100).';

alpha = eigs(double(A.' * A), 1) + 1;
theta = rand(100, 1);
lambda = 1;
for iter=1:5000
    theta = wthresh(theta + (1/alpha) * A.' * (y - A*theta), 's', lambda/(2*alpha));
end

stem(x)
figure
stem(y)
figure
stem(theta)
