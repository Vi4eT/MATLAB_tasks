%% Метод Крамера
matrix = [5.9507  1.3468 -1.7517  0.3284 -1.7145 -3.761  -4.2671 -0.1254 -1.8324 -3.2423
		 -2.4151 -0.0669  0.2975  1.2606  0.8317 -0.728   0.7213  2.0279  2.0937 -1.4506
		 -4.4816 -1.1502  9.0776  0.622  -0.4998 -3.0618 -1.5944  2.6214  2.2478  0.3083
		 -1.5108 -0.3819 -1.1529 -0.026   0.7621  0.9624 -0.5078  0.1086  2.2319  0.7678
		 -3.2314 -0.8235 -2.2375 -1.3416  5.5236 -0.5562 -1.5065  0.0295  0.9993  0.7226
		  1.092   0.3294 -1.0033  0.4390 -1.8597 10.2924 -4.6083 -9.3956 -2.3187  3.0877
		 -1.7932 -0.454  -1.3442 -0.7453 -1.7399 -1.3441 13.0325  0.5463  0.8250  0.8414
		 -0.8168 -0.2279  0.1156 -0.3336  0.3739 -4.2628  2.1442  5.0571  4.5388  2.1597
		 -0.0982 -0.0282  0.0411 -0.0399  0.0886 -0.6491  0.3831  0.7411 15.0656  2.5683
       	  0.3482  0.106  -0.3538  0.1397 -0.6473  3.3443 -2.3153 -4.141  -1.5785  0.9178];
rightpart = [12.1731, -10.9282, -6.934, 3.6014, 15.6151, 51.314, -72.5419, -28.1283, -46.1519, 23.9508];
A1 = matrix;
d = det(matrix);												%определитель
xmatrix = zeros(10);											%предварительное выделение памяти
for i = 1:10
	matrix(:, i) = rightpart;									%замена строки
	xmatrix(i) = det(matrix) / d;								%решение
	matrix = A1;												%возврат к первоначальному виду матрицы
end
xmatrix = xmatrix.';

%% Использование обратной матрицы
xmatrix1 = inv(matrix) * rightpart.';							%#ok<MINV>

%% Операция "\"
x = matrix \ rightpart.';

%% Сравнение решений
nev1 = matrix * xmatrix - rightpart.';							%вектор невязки, метод Крамера
nev2 = matrix * xmatrix1 - rightpart.';							%вектор невязки, использование обратной матрицы
nev3 = matrix * x - rightpart.';								%вектор невязки, операция "\"
m1 = norm(nev1);												%нормы
m2 = norm(nev2);
m3 = norm(nev3);

%% Погрешности
rightpart2 = rightpart + (rightpart .* rand(1, 10) * 0.02);		% 2%
x2 = matrix \ rightpart2.';										%решение
dy2 = norm(x - x2) / norm(x2);									%погрешность решения
dx2 = norm(rightpart2 - rightpart) / norm(rightpart);			%погрешность исходных данных
rightpart3 = rightpart + (rightpart .* rand(1, 10) * 0.03);		% 3%
x3 = matrix \ rightpart3.';										%решение
dy3 = norm(x - x3) / norm(x3);									%погрешность решения
dx3 = norm(rightpart3 - rightpart) / norm(rightpart);			%погрешность исходных данных
rightpart4 = rightpart + (rightpart .* rand(1, 10) * 0.04);		% 4%
x4 = matrix \ rightpart4.';										%решение
dy4 = norm(x - x4) / norm(x4);									%погрешность решения
dx4 = norm(rightpart4 - rightpart) / norm(rightpart);			%погрешность исходных данных
rightpart5 = rightpart + (rightpart .* rand(1, 10) * 0.05);		% 5%
x5 = matrix \ rightpart5.';										%решение
dy5 = norm(x - x5) / norm(x5);									%погрешность решения
dx5 = norm(rightpart5 - rightpart) / norm(rightpart);			%погрешность исходных данных
X = [dx2, dx3, dx4, dx5];
Y = [dy2, dy3, dy4, dy5];
disp('Значение cond')
c = cond(matrix);
disp(c)
Yc = Y * c;
digits(5);
Xv = vpa(X');
Yv = vpa(Y');
Ycv = vpa(Yc');
disp('Погр. данных Погр. решений Погр. данных*cond')
tab = [Xv, Yv, Ycv];
disp(tab)
plot(X, Y, 'r*')
title('График погрешностей');
xlabel('Погрешность данных');
ylabel('Погрешность решений');
grid on

%% Матрица СЛАУ
%% Полином
p = poly(matrix);												%формирование полинома
x = -2 : 0.01 : 16;												%аргументы полинома
y = polyval(p, x);												%значение полинома
figure
plot(x, y)														%построение графика
grid on
zero = roots(p);
ind = zero(imag(zero) == 0);									%отбор действительных корней
r = ind;
disp(r)
hold on
plot(r, 0, 'ro')												%построение корней
title('Полином')												%подпись графика
xlabel('x')														%подписи осей
ylabel('y')

%% Собственные числа, вектор невязки
[c, X1] = eig(matrix);											%собственные числа
for i = 1 : 10
	x = X1(i, i);
	c1 = c(:, i);
	d = matrix * c1 - x * c1;									%вектор невязки для каждой пары
end

%% Сравнение собственных чисел
rr = sort(zero);												%сортировка корней
ss = sort(diag(X1));											%сортировка значений
otn = max(abs((rr - ss) ./ rr));								%относительная погрешность