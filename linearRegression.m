%Paolo G. Peani and James W. Wedum

function [A, E] = linearRegression(T, K)
%A script to compute the linear regression line for the velocity data from
%class along with the correlation coefficient and a plot of the data.

%Start by performing some mathematical manipulations on T and K so that the
%resulting slope and y-intercept will correspond to the linearized version
%of the Arrhenius equation: ln(K) = -E/(R*T) + ln(A). Here ln(K)
%corresponds to y, -E/R is the slope, 1/T is x, and ln(A) is b
x = 1./T;
y = log(K);

%Also define R, the ideal gas constant:
R = 8.314;

%We'll need a few values in our computations.
xbar=mean(x);
ybar=mean(y);
n=size(x);
n=n(1);

%Compute the slope and y-intercept using the 
%formulas from class.
m=(n*sum(x.*y)-sum(x)*sum(y))/(n*sum(x.*x)-(sum(x))^2);
b=ybar-m*xbar;

%Now we need to extract the frequency factor, A
%and the activation energy, E
% b = ln(A) and m = -E/R so:
A = exp(b);
E = -m*R;
end

