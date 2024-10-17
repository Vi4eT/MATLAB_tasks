function y1 = fun1(x)
	a = 1.6;
	b = 0.5;
	c = 2.5;
	d = 0.4;
	y1 = exp(a + b * x .^ 2) .* sin(c + d * x);
end