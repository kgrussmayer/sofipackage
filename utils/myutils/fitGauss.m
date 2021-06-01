function  [res,f,x]= fitGauss(sig)

ft = fittype('a*exp(-((x-x0)^2/(2*s^2)))+b');
[A,pos] = max(sig);
b = min(sig);
t = sig-b;t = t./max(t);
s = sum(t > 0.5);
[f,x] = fit((1:length(sig))',sig',ft,'Startpoint',[A b s pos]);

res = sqrt(8*log(2))*f.s;