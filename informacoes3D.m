function informacoes3D(Lx,Ly,Lz,nx,ny,nz,v,d,BETA,ETA,nt)
disp('------------------------------')
disp('----- Covariance function ----')
if nt==1
    disp('Exponential fields')
    disp('Correlation length (x):')
ETA
end
if nt==3
    disp('Square exponential fields')
    disp('Correlation length (x):')
ETA
end
if nt==2
    disp('Fractal fields')
  disp('Hurst coefficient:')
BETA
end
disp('------------------------------')
disp('-------- DOMAIN SIZE ---------')
Lx
Ly
Lz
disp('------------------------------')
disp('------------ MESH ------------')
nx
ny
nz
%
n=length(d);
if(n==0)
disp('------------------------------')
disp('------ NO CONDITIONING -------') 
else
disp('------------------------------')
disp('-------- CONDITIONING --------')
disp('---- Number of points (n) ----')
  n
disp('------ Positions (x,y) -------')
  v
disp('---------- Values ------------')
  d
disp('------------------------------')
end
disp('##############################')
disp('')
