function y = pol(x)
	a = 0.7;
	b = 3.43;
	c = 5.516;
	d = 2.912;
	p = [a, b, c, d];
	y = polyval(p, x);
end