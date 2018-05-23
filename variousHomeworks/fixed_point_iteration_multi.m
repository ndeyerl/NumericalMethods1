function r=fixed_point_iteration_multi(x0 , y0, precision)
xi=x0;
yi=y0;
xii=(8/yi)-sin(xi)/(yi.^2);
yii=3*xi-(1/2)*cos(xi);

while abs(xii-xi)>precision || abs(yii-yi) > precision %absolute precision
    xi=xii;
    yi=yii;
    xii=xii(xi,yi);
    yii=yii(xi,yi);
end
r= [xii, yii]
end