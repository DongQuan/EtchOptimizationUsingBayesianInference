function F = CalculateRates(A,B,Te)
k = zeros(length(A),1);

for i=1:length(A)
    k(i) = A(i)*exp(-B(i)/(Te)); %Arrhenius equation
end

%make dimensionless
F = k;