% Paolo G. Peani and James W. Wedum

function [A,E,b] = multipleLinearRegression(T, k)
%Uses method of multiple linear regression to find linear
% fit of Arrhenius Equation. k=A*(T^b)*(e^(-E/(RT)))
%   @x = x coordinates of data points
%   @y = y coordinates of data points
%
%   @A,E,b = constants for Arrhenius equation labled 
%       to match equation

%manipulate x and y into linear form:
%   ln(k) + e = ln(A) + b*ln(T) + E*(-1/RT) + e
e = exp(1);
y = log(k);
x1 = log(T);
x2 = -1./(8.314.*T);

%find matrix values
%|  n    sumx1   sumx2 |   |a0|   | sumy |
%|sumx1 sumx1Sq sumx1x2| * |a1| = |sumx1y|
%|sumx2 sumx1x2 sumx2Sq|   |a2|   |sumx2y|
n=length(k);
sumx1=sum(x1);
sumx2=sum(x2);
sumx1Sq = sum(x1.^2);
sumx2Sq = sum(x2.^2);
sumx1x2=sum(x1.*x2);
sumy=sum(y);
sumx1y=sum(x1.*y);
sumx2y=sum(x2.*y);

%set up and solve system of linear equations
%innitial guesses of 1 for a0, 0 for a1 and a2
a0=1;
a1=0;
a2=0;
loopcount = 0;

%very close approximation of solutions, though not likely exact
while((n*a0+sumx1*a1+sumx2*a2~=sumy ...
    && sumx1*a0+sumx1Sq*a1+sumx1x2*a2~=sumx1y ...
    && sumx2*a0+sumx1x2*a1+sumx2Sq*a2~=sumx2y)...
    ||loopcount<10000)
a0=(sumy-sumx1*a1-sumx2*a2)/n;
a1=(sumx1y-sumx1*a0-sumx1x2*a2)/sumx1Sq;
a2=(sumx2y-sumx2*a0-sumx1x2*a1)/sumx2Sq;

loopcount=loopcount+1;
end


A=exp(a0);
E=a2;
b=a1;

%Display graph and data for testing, uncomment for graph
%----------------------------------
% hold('on')
% plot(T,k,'*g')
% plot(T,A.*T.^b.*exp(-E./(8.314.*T)))
% hold('off')
%----------------------------------

%Find and display rms value
St=sum((k-mean(k)).^2);
S=sum((k-(A.*T.^b.*exp(-E./(8.314.*T)))).^2);
r=sqrt((St-S)/St);
sprintf('RMS deviation for multiple linear regression model: %.3f',r)

end

