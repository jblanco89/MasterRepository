format long;
clear, clc
% Problem:
% Calculate revolution solid volumen generated spining around axis OX of
% this function. Use Trapezium, midpoint, simpson 1/3 and simpson 1/8 with
% h = 0.1. Problem taken from Zaragoza University. 

%---------------%
%| y = 1 - x^2;|%
%| x = -1;     |%
%| x = 1;      |%
%---------------%


h = 0.1;

f = @(x) 1 - x.^2;
x1a = -1;
x1b = 1;
x = linspace(x1a, x1b, 20);
y = f(x);

figure;
plot(x, y, 'LineWidth', 2, 'DisplayName', '1 - x^2');
hold on;
xline(x1a,'r','LineWidth', 1.5,'Label', '-1');
xline(x1b,'r','LineWidth', 1.5,'Label', '1');
xlim([-1.20 1.20])
fill([x, fliplr(x)], [y, zeros(size(y))], 'c', ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none', 'DisplayName', ...
    'Área del sólido acotado');
grid on;
legend('show');

% Solid Volume is:
% V = pi * I(f(x)^2, a, b)dx
f2 = @(x) (1 - x.^2).^2;
y2 = f2(x);

n = (x1b - x1a) / h;

% MidPoint

[Imid, Emid] = midpoint(f2, 0, x1a, x1b, n);
Vmid = pi * Imid;

% Trapezium
[Itra, Etra] = trapezium(f2, 0, x1a, x1b, n);
Vtra = pi* Itra;

% Simpson 1/3
[Isimp13, Esimp13] = simpson13(f2, 0, x1a, x1b, n);
% [Isimp13, Esimp13] = simpson13(0, y2, x1a, x1b, n);
Vsimp13 = pi* Isimp13;

% Simpson 3/8
[Isimp38, Esimp38] = simpson38(f2, 0, x1a, x1b, n-2);
% [Isimp38, Esimp38] = simpson38(0, y2(1:end-1), x1a, x1b, n-2);
Vsimp38 = pi* Isimp38;

Ints = [Vtra, Vsimp13, Vsimp38, Vmid];
Errs = [Etra, Esimp13, Esimp38, Emid];

result = table(Ints', Errs','RowNames', ...
    {'Trapezium', 'Simpson 1/3', 'Simpson 3/8', 'MidPoint'}, ...
    'VariableNames',{'Volume', 'Error'});

disp(result);

