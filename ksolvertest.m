function F = ksolvertest(x,Act,B)
scale = 1;
F = [
x(1)/scale - Act(1)*exp(-B(1)/x(8));%k1
x(2)/scale - Act(2)*exp(-B(2)/x(8));%k2
x(3)/scale - Act(3)*exp(-B(3)/x(8));%k3
x(4)/scale - Act(4)*exp(-B(4)/x(8));%k4
x(5)/scale - Act(5)*exp(-B(5)/x(8));%k5
x(6)/scale - Act(6)*exp(-B(6)/x(8));%k6
x(7)/scale - Act(7)*exp(-B(7)/x(8))];
end