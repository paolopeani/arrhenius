% MAIN SCRIPT
% Paolo G. Peani and James W. Wedum

%Start by initializing known values:
%R, the ideal gas constant, is pre-defined:
R=8.314;

%T_a, the temperature, and K, the reaction rate are given as 
%experimental data, so we need to load them in:
data=load('data.txt');
T_a = data(:,1);
K = data(:,2);

%We will also want the mean of K for use in calculating the correlation
%coefficient:
Kbar = mean(K);

%We now pass T_a and K to our subfunctions.
%group2LinearRegression will return the A and E values for the simplified
%Arrhenius equation. group2MultipleLinearRegression will return the A, E,
%and b values for the more sophisticated Arrhenius equation.
[A1,E1]=linearRegression(T_a,K);
[A2,E2,b2]=multipleLinearRegression(T_a,K);

%Once we have found the unknown parameters, we can enter them into our
%equations to determine the K values for our models. K1 is a vector
%representing the K values for the simplified Arrhenius fit. K2 is a vector
%representing the K values for the sophisticated Arrhenius fit.
K1 = A1.*exp((-E1)./(R.*T_a));
K2 = A2.*(T_a.^b2).*exp(-E2./(R.*T_a));

%Now calculate and display the correlation coefficients for both models.
%First model:
r1=sqrt((sum((K-Kbar).^2)-sum((K-K1).^2))/(sum((K-Kbar).^2)));
str1 = sprintf('The correlation coefficient for the simplified Arrhenius Equation is: %f', r1);
disp(str1);

%Second Model:
r2=sqrt((sum((K-Kbar).^2)-sum((K-K2).^2))/(sum((K-Kbar).^2)));
str2 = sprintf('The correlation coefficient for the sophisticated Arrhenius Equation is: %f', r2);
disp(str2);

%Now plot our fits against the experimental data:
subplot(2,1,1);
plot(T_a,K,'.',T_a,K1);
title('Linear Fit for Simplified Arrhenius Equation');
xlabel('Temperature (Kelvin)');
ylabel('Reaction Rate');
legend('Data', 'Fit');


subplot(2,1,2);
plot(T_a,K,'.',T_a,K2);
title('Multiple Linear Fit for Sophisticated Arrhenius Equation');
xlabel('Temperature (Kelvin)');
ylabel('Reaction Rate');
legend('Data', 'Fit');


