

fy(x)=P0y * exp(x * P1y) ;
P0y=1
P1y=0.001
fit fy(x) 'test_non_linear_fit_sigmax_gnuplot.data' i 0 u   1:2:4 yerrors via P0y,P1y


f(x)=P0 * exp(x * P1) ;
P0=1
P1=0.001
fit f(x) 'test_non_linear_fit_sigmax_gnuplot.data' i 0 u   1:2:3:4 xyerrors via P0,P1


fc(x)=1.10204e-05 * exp(x * 0.0427657) ;
fcy(x)=1.75235e-05 * exp(x * 0.0414713) ;


p [305:] 'test_non_linear_fit_sigmax_gnuplot.data' u 1:2:3:4 ps 1.3 pt 2 t 'Data' w xyerrorbars,\
f(x) t 'Fit, xy error',   fy(x) t 'Fit, y error' ,\
fc(x) t 'Fit, C++ xy error ' , fcy(x) t 'Fit, C++ y error '
