function z = Cov3Df(xv,yv,sigma,eta1,eta2,eta3)
% x = (xv(:,1)-yv(:,1));
% y = (xv(:,2)-yv(:,2));
% funcao de covariancia 2D 
% Note que xv = (x1,x2), yv = (y1,y2)
z = sigma*exp(-((xv(:,1)-yv(:,1)).^2)/(2*eta1) -...
    ((xv(:,2)-yv(:,2)).^2)/(2*eta2)-((xv(:,3)-yv(:,3)).^2)/(2*eta3));
%z = sigma*exp(-((x).^2)/(2*eta1) - ((y).^2)/(2*eta2));