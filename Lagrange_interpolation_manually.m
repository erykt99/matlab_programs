X = [10, 20, 30, 40, 50];
Y = [0.98, 0.93, 0.86, 0.76, 0.64];
figure(1)
plot(X, Y);
hold on

x0 = 21;

t = 0;
y  = 0.0;

for k = 1:1:length(X)
t = 1;
for j = 1:1:length(X)
if(j ~= k)
    t = t *((x0 - X(j))/(X(j) - X(k)));
end
end
    y = y + (t*Y(k)); 
end

plot(x0,y,'+');
title("Lagrange")
figure(2)
x1 = 1:1:50;
y1 = interp1(X, Y, x1);
plot(x1, y1);
title("Matlab fun")

