% x = [-3:.1:3];
% N=1;
% norm = N*normpdf(x,0,1);
% figure;
% plot(norm,x)

r1 = rand(1000,1);
f = histogram(r1);%# create histogram from a normal distribution.
f.Normalization = 'probability';
