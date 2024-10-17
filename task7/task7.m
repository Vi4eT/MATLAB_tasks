%% Пределы
%% 428
syms x m n
disp('№428:')
s428 = (m / (1 - x ^ m)) - (n / (1 - x ^ n));
pretty(s428)														%красивый вывод
s428lim = limit(s428, x, 1);										%вычисление предела
disp('№428 ответ:')
pretty(s428lim)														%красивый вывод

%% 488
syms a
disp('№488:')
s488 = (sin(a + 2 * x) - 2 * sin(a + x) + sin(a)) / x ^ 2;
pretty(s488)														%красивый вывод
s488lim = limit(s488, x, 0);										%вычисление предела
disp('№488 ответ:')
pretty(s488lim)														%красивый вывод

%% 596(б)
disp('№596(б):')
s596 = 1 / (1 + exp(1 / x));
pretty(s596)														%красивый вывод
s596lim = limit(s596, x, 0, 'right');								%вычисление предела
disp('№596(б) ответ:')
pretty(s596lim)														%красивый вывод

%% Производные
%% 892
disp('№892:')
s892 = (1 / (2 * sqrt(6))) * log((x * sqrt(3) - sqrt(2)) / (x * sqrt(3) + sqrt(2)));
pretty(s892)														%красивый вывод
s892dif = diff(s892, 'x', 1);										%вычисление произодной
disp('№892 ответ:')
pretty(simplify(s892dif))											%красивый вывод упрощенного ответа

%% 919
disp('№919:')
s919 = x * asin(sqrt(x / (1 + x))) + atan(sqrt(x)) - sqrt(x);
pretty(s919)														%красивый вывод
s919dif = diff(s919, 'x', 1);										%вычисление произодной
disp('№919 ответ:')
pretty(simplify(s919dif))											%красивый вывод

%% Предел
%% 1335
disp('№1335:')
s1335 = (log(sinh(x) + sqrt(1 + sinh(x) ^ 2)) - log(sin(x) + sqrt(1 + sin(x) ^ 2))) / (sinh(x) - sin(x));
pretty(s1335)														%красивый вывод
s1335lim = limit(s1335, x, 0);										%вычисление предела
disp('№1335 ответ:')
pretty(s1335lim)													%красивый вывод

%% Кубическое уравнение
p = [1, -3, 6, -2];													%коэффициенты полинома
X1 = roots(p);														%корни полинома (roots)
X1real = X1(imag(X1) == 0);											%вещественные корни полинома (roots)
f = str2sym('x ^ 3 - 3 * x ^ 2 + 6 * x - 2');
X2 = solve(f);														%аналитическое решение уравнения
X2real = double(X2);
X2rreal = X2real(imag(X2real) == 0);								%вещественные корни полинома
k1 = vpa(X1real);													%roots
k2 = vpa(X2rreal);													%solve

%% Ряд Тейлора
syms x;
Y = 2 * x + log10(x) + 0.5;
t2 = taylor(Y, x, 2, 'order', 2);									%разложение в ряд Тейлора, 2 слагаемых
t3 = taylor(Y, x, 2, 'order', 3);									%3
t4 = taylor(Y, x, 2, 'order', 4);									%4
t5 = taylor(Y, x, 2, 'order', 5);									%5
pt2 = sym2poly(t2);													%коэффициенты полинома разложения
pt3 = sym2poly(t3);
pt4 = sym2poly(t4);
pt5 = sym2poly(t5);
X = linspace(0.1, 3.9, 100);										%интервал построения графика
y2 = polyval(pt2, X);												%значения полиномов
y3 = polyval(pt3, X);
y4 = polyval(pt4, X);
y5 = polyval(pt5, X);
subplot(2, 2, 1)
Y1 = subs(Y, x, X);
plot(X, Y1)
hold on
plot(X, y2, 'r')													%построение полинома из 2 слагаемых
grid on
title('Ряд Тейлора (2 слагаемых)');
xlabel('x')
ylabel('y')
subplot(2, 2, 2)
plot(X, Y1)
hold on
plot(X, y3, 'r')													%построение полинома из 3 слагаемых
grid on
title('Ряд Тейлора (3 слагаемых)');
xlabel('x')
ylabel('y')
subplot(2, 2, 3)
plot(X, Y1)
hold on
plot(X, y4, 'r')													%построение полинома из 4 слагаемых
grid on
title('Ряд Тейлора (4 слагаемых)');
xlabel('x')
ylabel('y')
subplot(2, 2, 4)
plot(X, Y1)
hold on
plot(X, y5, 'r')													%построение полинома из 5 слагаемых
grid on
title('Ряд Тейлора (5 слагаемых)');
xlabel('x')
ylabel('y')

%% Решение СЛАУ
a = rand(3);														%матрица 3 порядка из случайных чисел
A = round(-9 + a * 18);												%приведение к интервалу [-9;9]
y = rand(3, 1);
Y = round(-9 + y * 18);												%правая часть матрицы
while (det(A) == 0)													%проверка существования и единственности решения 
	A(1) = -9 + 18 * rand();
end
syms x1 x2 x3
y1 = A(1, 1) * x1 + A(1, 2) * x2 + A(1, 3) * x3 - Y(1);
y2 = A(2, 1) * x1 + A(2, 2) * x2 + A(2, 3) * x3 - Y(2);				%СЛАУ
y3 = A(3, 1) * x1 + A(3, 2) * x2 + A(3, 3) * x3 - Y(3);
disp(y1);
disp(y2);
disp(y3);
X = inv(A) * Y;														%#ok<MINV> %решение с помощью обратной матрицы
digits(10);
X = vpa(X);
rot = solve(y1, y2, y3);											%аналитическое решение
XX(1, 1) = vpa(rot.x1);
XX(2, 1) = vpa(rot.x2);
XX(3, 1) = vpa(rot.x3);
del = abs(XX-X);													%разница измерений
rez = [XX, X, del];
disp('  Обр. матрица | Аналитически | Разница')
disp(rez)