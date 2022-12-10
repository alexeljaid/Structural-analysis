% Linear regression script for the calculation of alpha and beta

clear all
close all

Es=2000; % MPa
rhobar=[0.05 0.1 0.15 0.2 0.25]; % relative density

% Ex=[1. 2. 3. 4. 5.]; %This is going to be our different values of Ex and Ey
Ex=[0.3750 2.9999 10.1245 23.9978 46.8684]; %Ex for hexagon
% Ex=[33.3356 66.6852 100.0625 133.4814 166.9558]; %Ex for triangle (Option 1)
% Ex=[33.3333 66.6667 100.0000 133.3333 166.6667]; %Ex for triangle (Option 2)
% Ey=[0.3750 3.0000 10.1248 23.9993 46.8728]; %Ey for hexagon
% Ey=[33.3356 66.6852 100.0625 133.4814 166.9558]; %Ey for triangle (Option 1)
% Ey=[33.3339 66.6713 100.0156 133.3703 166.7389]; %Ey for triangle (Option 2)

for n=1:5
    logrhobar(n)=log(rhobar(n));
    logEx(n)=log(Ex(n));
end

xabsis=linspace(min(logrhobar),0,10);


p = polyfit(logrhobar,logEx,1);
lininterp = polyval(p,xabsis);

a = lininterp(length(xabsis));
disp(a)
alpha = exp(a-log(Es));


beta = (lininterp(2)-lininterp(1))/(xabsis(2)-xabsis(1));


plot(xabsis,lininterp)
hold on
scatter(logrhobar,logEx)
xlabel("log(rhobar)")
ylabel("log(Ex)")
title("Alpha="+num2str(alpha)+", Beta="+num2str(beta))
hold off
