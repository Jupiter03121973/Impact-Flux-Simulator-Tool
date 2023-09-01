% Gegebene Wertepaare
x = [6975, 6875, 6775,6675];
y1 = [0.010543478, 0.016640895, 0.031439153, 0.045951505];
y2 = [0.008531232, 0.014430455, 0.014423032, 0.03229449];
y3 = [0.007947652, 0.007794118, 0.014833843, 0.017112035];
y4 = [0.004997146, 0.006376737, 0.010627045, 0.007517518];
y5 = [0.003356826, 0.004332119, 0.008921265, 0.023702337];
y6 = [0.00722997,  0.007227208, 0.008572207, 0.011202849];

% Polynom 3. Grades anpassen
p1 = polyfit(x, y1, 3);
p2 = polyfit(x, y2, 3);
p3 = polyfit(x, y3, 3);
p4 = polyfit(x, y4, 3);
p5 = polyfit(x, y5, 3);
p6 = polyfit(x, y6, 3);

matrix = vertcat(p1, p2, p3, p4, p5, p6);
% Matrix in Exponential-Schreibweise ausgeben
for i = 1:size(matrix, 1)
    fprintf('%1.3e,  %1.3e;\n', matrix(i, 1), matrix(i, 2),matrix(i, 3),matrix(i, 4));
end
% Wert fuer SMA
SMA = 6800;

% Ausgleichskurve berechnen
f = p1(1) * SMA.^3 + p1(2) * SMA.^2 + p1(3) * SMA + p1(4)

% Ausgleichskurve berechnen
f = polyval(p1, SMA);

disp(['f(SMA) = ', num2str(f)]);

% Ausgabe des Polynoms mit eingesetzten Koeffizienten
disp(['f(SMA) = ', num2str(p1(1)), ' * SMA.^3 + ', num2str(p1(2)), ' * SMA.^2 + ', num2str(p1(3)), ' * SMA + ', num2str(p1(4))]);