%% ������������ ��������
a = 0.7;
b = 3.43;
c = 5.516;
d = 2.912;
p = [a, b, c, d];							%���� ������������� ��������
x = linspace(-2.1, -1.18, 100);				%��������� ��������
y = polyval(p, x);							%�������� ��������
r = roots(p);								%����� ��������
x0 = r(imag(r) == 0);						%�������������� ����� ��������
p1 = polyder(p);							%����������� ��������
r1 = roots(p1);								%����� �����������
x2 = r1(imag(r1) == 0);						%�������������� ����� �����������
y2 = polyval(p, x2);						%�������� �����������
p2 = polyder(p1);							%����������� ������� �������
n = length(p2);								%���-�� ������������� ����������� ������� �������

%% ������������ �������
x1 = linspace(-1, 2, 100);					%��������� �������
y1 = fun1(x1);								%�������� �������

%% ������ 1.1
subplot(1, 2, 1);							%���������� ��������
plot(x, y);
hold on
plot(x0, 0, 'rs')							%���������� ������ ��������
for i = 1 : n
	p3 = polyval(p2, x2(i));				%�������� ����������� ������� �������
	hold on
	if (p3 < 0)
		plot(x2(i), y2(i), 'ro');			%���������� ���������
		text(-1.86, 0.016, 'max');			%������� ���������
	else
		plot(x2(i), y2(i), 'ro');			%���������� ��������
		text(-1.45, -0.01, 'min');			%������� ��������
	end
end
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('�������');							%������� �������
text(-2.02, 0.0015, 'root1');				%������� ������
text(-1.62, 0.0015, 'root2');
text(-1.32, 0.0015, 'root3');

%% ������ 1.2
subplot(1, 2, 2);							%���������� �������
plot(x1, y1);
hold on
[m1, n1] = fzero(@fun1, [1, 2]);			%������ �������
plot(m1, n1, 'or');							%���������� ����� �������
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('�������');							%������� �������
text(1.545, 0.35, 'root');					%������� �����

%% ���� 2
%% ������ 2.1
figure										%������� � ������ ����
subplot(2, 2, 1);							%���������� �������� 1
plot(x, y);
hold on
X = linspace(-2.1, -1.18, 7);
Y = polyval(p, X);
plot(X, Y, 'ro')
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 1');							%������� �������

%% ������ 2.2
subplot(2, 2, 2);							%���������� �������� 2
plot(x, y, 'r:');
hold on
plot(X, Y, 'bx')
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 2');							%������� �������

%% ������ 2.3
subplot(2, 2, 3);							%���������� �������� 3
plot(x, y, 'k--');
hold on
plot(X, Y, 'g*')
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 3');							%������� �������

%% ������ 2.4
subplot(2, 2, 4);							%���������� �������� 4
plot(x, y, 'g-.');
hold on
plot(X, Y, 'kh')							%���������� ������ �������� 4
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 4');							%������� �������

%% ���� 3
%% ������ 3.1
figure										%������� � ������ ����
subplot(2, 2, 1);							%���������� ������� 1
plot(x1, y1);
hold on
X1 = linspace(-1, 2, 7);
Y1 = fun1(X1);
plot(X1, Y1, 'or');
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 1');							%������� �������

%% ������ 3.2
subplot(2, 2, 2);							%���������� ������� 2
plot(x1, y1, 'r:');
hold on
plot(X1, Y1, 'bx');							%���������� ����� ������� 2
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 2');							%������� �������

%% ������ 3.3
subplot(2, 2, 3);							%���������� ������� 3
plot(x1, y1, 'k--');
hold on
plot(X1, Y1, 'g*');
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 3');							%������� �������

%% ������ 3.4
subplot(2, 2, 4);							%���������� ������� 4
plot(x1, y1, 'g-.');
hold on
plot(X1, Y1, 'kh');							%���������� ����� ������� 4
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 4');							%������� �������

%% ���� 4
%% ������ 4.1
figure										%������� � ������ ����
subplot(2, 2, 1);							%���������� �������� 1
x = linspace(-2.5, 0.5, 100);				%��������� �������� ������� ��������
y = polyval(p, x);							%�������� ��������
plot(x, y);									%���������� �������� 1
hold on
plot(x0, 0, 'rs')							%���������� ������ �������� 1
for i = 1:n
	p3 = polyval(p2, x2(i));				%�������� ����������� ������� �������
	hold on
	if (p3 < 0)
		plot(x2(i), y2(i), 'ro');			%���������� �����������
	end
end
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 1');							%������� �������
hold on
plot(x1, y1);								%���������� ������� 1
hold on
plot(m1, n1, 'ro');							%���������� ����� ������� 1

%% ������ 4.2
subplot(2, 2, 2);							%���������� �������� 2
plot(x, y, 'r:');							%���������� �������� 2
hold on
plot(x0, 0, 'bx')							%���������� ������ �������� 2
for i = 1:n
	p3 = polyval(p2, x2(i));				%�������� ����������� ������� �������
	hold on
	if (p3 < 0)
		plot(x2(i), y2(i), 'b.');			%���������� �����������
	end
end
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 2');							%������� �������
hold on
plot(x1, y1, 'r:');							%���������� ������� 2
hold on
plot(m1, n1, 'bx');							%���������� ����� ������� 2

%% ������ 4.3
subplot(2, 2, 3);							%���������� �������� 3
plot(x, y, 'k--');							%���������� �������� 3
hold on
plot(x0, 0, 'g*')							%���������� ������ �������� 3
for i = 1:n
	p3 = polyval(p2, x2(i));				%�������� ����������� ������� �������
	hold on
	if (p3 < 0)
		plot(x2(i), y2(i), 'gd');			%���������� �����������
	end
end
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 3');							%������� �������
hold on
plot(x1, y1, 'k--');						%���������� ������� 3
hold on
plot(m1, n1, 'g*');							%���������� ����� ������� 3

%% ������ 4.4
subplot(2, 2, 4);							%���������� �������� 4
plot(x, y, 'g-.');							%���������� �������� 4
hold on
plot(x0, 0, 'kh')							%���������� ������ �������� 4
for i = 1:n
	p3 = polyval(p2, x2(i));				%�������� ����������� ������� �������
	hold on
	if (p3 < 0)
		plot(x2(i), y2(i), 'rp');			%���������� �����������
	end
end
grid on										%�����
xlabel('x');								%������� ����
ylabel('y');
title('������� 4');							%������� �������
hold on
plot(x1, y1, 'g-.');						%���������� ������� 4
hold on
plot(m1, n1, 'kh');							%���������� ����� ������� 4