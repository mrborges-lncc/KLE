function cc = Cov3D(xv,yv,sigma,eta1,eta2,eta3)
 x = (xv(:,1)-yv(:,1));
 y = (xv(:,2)-yv(:,2));
 z = (xv(:,3)-yv(:,3));
% funcao de covariancia 2D 
% Note que xv = (x1,x2), yv = (y1,y2)
% z = sigma*exp(-((xv(:,1)-yv(:,1)).^2)/(eta1*eta1) - ((xv(:,2)-yv(:,2)).^2)/(eta2*eta2));
cc = sigma*exp(-sqrt((x.^2)/(eta1) + (y.^2)/(eta2)+ (z.^2)/(eta3)));