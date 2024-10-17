%% �������
%% 428
syms x m n
disp('�428:')
s428 = (m / (1 - x ^ m)) - (n / (1 - x ^ n));
pretty(s428)														%�������� �����
s428lim = limit(s428, x, 1);										%���������� �������
disp('�428 �����:')
pretty(s428lim)														%�������� �����

%% 488
syms a
disp('�488:')
s488 = (sin(a + 2 * x) - 2 * sin(a + x) + sin(a)) / x ^ 2;
pretty(s488)														%�������� �����
s488lim = limit(s488, x, 0);										%���������� �������
disp('�488 �����:')
pretty(s488lim)														%�������� �����

%% 596(�)
disp('�596(�):')
s596 = 1 / (1 + exp(1 / x));
pretty(s596)														%�������� �����
s596lim = limit(s596, x, 0, 'right');								%���������� �������
disp('�596(�) �����:')
pretty(s596lim)														%�������� �����

%% �����������
%% 892
disp('�892:')
s892 = (1 / (2 * sqrt(6))) * log((x * sqrt(3) - sqrt(2)) / (x * sqrt(3) + sqrt(2)));
pretty(s892)														%�������� �����
s892dif = diff(s892, 'x', 1);										%���������� ����������
disp('�892 �����:')
pretty(simplify(s892dif))											%�������� ����� ����������� ������

%% 919
disp('�919:')
s919 = x * asin(sqrt(x / (1 + x))) + atan(sqrt(x)) - sqrt(x);
pretty(s919)														%�������� �����
s919dif = diff(s919, 'x', 1);										%���������� ����������
disp('�919 �����:')
pretty(simplify(s919dif))											%�������� �����

%% ������
%% 1335
disp('�1335:')
s1335 = (log(sinh(x) + sqrt(1 + sinh(x) ^ 2)) - log(sin(x) + sqrt(1 + sin(x) ^ 2))) / (sinh(x) - sin(x));
pretty(s1335)														%�������� �����
s1335lim = limit(s1335, x, 0);										%���������� �������
disp('�1335 �����:')
pretty(s1335lim)													%�������� �����

%% ���������� ���������
p = [1, -3, 6, -2];													%������������ ��������
X1 = roots(p);														%����� �������� (roots)
X1real = X1(imag(X1) == 0);											%������������ ����� �������� (roots)
f = str2sym('x ^ 3 - 3 * x ^ 2 + 6 * x - 2');
X2 = solve(f);														%������������� ������� ���������
X2real = double(X2);
X2rreal = X2real(imag(X2real) == 0);								%������������ ����� ��������
k1 = vpa(X1real);													%roots
k2 = vpa(X2rreal);													%solve

%% ��� �������
syms x;
Y = 2 * x + log10(x) + 0.5;
t2 = taylor(Y, x, 2, 'order', 2);									%���������� � ��� �������, 2 ���������
t3 = taylor(Y, x, 2, 'order', 3);									%3
t4 = taylor(Y, x, 2, 'order', 4);									%4
t5 = taylor(Y, x, 2, 'order', 5);									%5
pt2 = sym2poly(t2);													%������������ �������� ����������
pt3 = sym2poly(t3);
pt4 = sym2poly(t4);
pt5 = sym2poly(t5);
X = linspace(0.1, 3.9, 100);										%�������� ���������� �������
y2 = polyval(pt2, X);												%�������� ���������
y3 = polyval(pt3, X);
y4 = polyval(pt4, X);
y5 = polyval(pt5, X);
subplot(2, 2, 1)
Y1 = subs(Y, x, X);
plot(X, Y1)
hold on
plot(X, y2, 'r')													%���������� �������� �� 2 ���������
grid on
title('��� ������� (2 ���������)');
xlabel('x')
ylabel('y')
subplot(2, 2, 2)
plot(X, Y1)
hold on
plot(X, y3, 'r')													%���������� �������� �� 3 ���������
grid on
title('��� ������� (3 ���������)');
xlabel('x')
ylabel('y')
subplot(2, 2, 3)
plot(X, Y1)
hold on
plot(X, y4, 'r')													%���������� �������� �� 4 ���������
grid on
title('��� ������� (4 ���������)');
xlabel('x')
ylabel('y')
subplot(2, 2, 4)
plot(X, Y1)
hold on
plot(X, y5, 'r')													%���������� �������� �� 5 ���������
grid on
title('��� ������� (5 ���������)');
xlabel('x')
ylabel('y')

%% ������� ����
a = rand(3);														%������� 3 ������� �� ��������� �����
A = round(-9 + a * 18);												%���������� � ��������� [-9;9]
y = rand(3, 1);
Y = round(-9 + y * 18);												%������ ����� �������
while (det(A) == 0)													%�������� ������������� � �������������� ������� 
	A(1) = -9 + 18 * rand();
end
syms x1 x2 x3
y1 = A(1, 1) * x1 + A(1, 2) * x2 + A(1, 3) * x3 - Y(1);
y2 = A(2, 1) * x1 + A(2, 2) * x2 + A(2, 3) * x3 - Y(2);				%����
y3 = A(3, 1) * x1 + A(3, 2) * x2 + A(3, 3) * x3 - Y(3);
disp(y1);
disp(y2);
disp(y3);
X = inv(A) * Y;														%#ok<MINV> %������� � ������� �������� �������
digits(10);
X = vpa(X);
rot = solve(y1, y2, y3);											%������������� �������
XX(1, 1) = vpa(rot.x1);
XX(2, 1) = vpa(rot.x2);
XX(3, 1) = vpa(rot.x3);
del = abs(XX-X);													%������� ���������
rez = [XX, X, del];
disp('  ���. ������� | ������������ | �������')
disp(rez)