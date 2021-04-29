clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    3D STOCHASTIC FIELDS GENERATOR       
%                   BASED ON KARHUNEN-LOEVE EXPANSION
% AUTHOR....: MARCIO RENTES BORGES
% EMAIL.....: marcio.rentes.borges@gmail.com
% DATE......: 11/11/2013
% REVIEWED..: 06/10/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputbox = 10; % if == 1 display a dialog box to input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart = tic;
home  ='/home/mrborges/MCMC_par/trunk/forecast/simuladorRigido/exp/fields/' % folder to save fields
home  ='/home/mrborges/Upscaling/KLgenerator/campos/' % folder to save fields
home  = '/home/mrborges/Dropbox/DeepLearning/Proj_simul/simuladorRigido/expUps_fine_scale/fields/'
home  ='/home/mrborges/MCMC_par/trunk/forecast/simuladorRigido/exp/fields/' % folder to save fields
home  = '~/MCMC_par/trunk/simuladorRigido/exp/fields/';
homeT = '~/MCMC_par/trunk/gera_KL/MATLAB/';
%home  ='/home/mrborges/Dropbox/ADWR/Estudos/fields/' % folder to save fields
%home  = '/home/mrborges/Dropbox/KLE/fields/';
homep ='/home/mrborges/Dropbox/ADWR/Estudos/paraview/' % folder to save fields
homep ='~/MCMC_par/trunk/simuladorRigido/exp/fields/'; % folder to save fields
home_fig= './pictures/' % folder to save pictures
home_fig = '/home/mrborges/Dropbox/DeepLearning/Proj_simul/simuladorRigido/expUps_fine_scale/paraview/'
home_fig = '~/MCMC_par/trunk/figuras/'
%home_fig = '/home/mrborges/Dropbox/ADWR/Estudos/figuras/'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntipo=3; 
%% define the covariance function: if ntipo == 1 exponential;
%%                                          == 2 fractal;
%%                                          == 3 square exponential
%%                                          == 4 Matern variogram
nu = 0.5; % Matern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
graf=1;  % if == 1 print the eigenvalues picture
fig =1;  % if == 1 print 3D field (last one)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% physical dimensions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx  = 1.00;
Ly  = 1.00;
Lz  = 1.00;
%%%%% mesh for covariance matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx  = 100;
ny  = 100;
nz  = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mesh for interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NX  = 40;
NY  = 40;
NZ  = 20;
interpolacao = 10; % if == 1 the eigenvector are interpolated to this mesh%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta1  = 0.20; % correlation length in the x direction
eta2  = 0.20; % correlation length in the y direction
eta3  = 0.02; % correlation length in the z direction
Nrand = 1;  % total number of realizations
M     = 0;  % number of terms used in the KL expansion. OBS: if == 0 it 
            % uses the maximum number of terms (nx^2 x ny^2 x nz^2)
TIPOINPUT = 10; % if == 1 reads the conditioned points from the file
                % indicated in "file_input_cond"
file_input_cond = '/home/mrborges/MCMC_par/trunk/gera_KL/FORTRAN_RW/in/input_cond.dat';
file_input_cond = '/home/mrborges/Dropbox/ADWR/Estudos/input_cond.dat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varY  = 1.0;           % field variance
beta  = 0.5;           % if ntipo == 2 this is the Hurst coefficient
cutoff= Lx/double(nx); % cutoff used when ntipo == 2
alpha = 1.0;           % KEEP == 1
tipo_prt  = 1;         % if == 1 print the fields in the LNCC format,
                       % if == 0 print in the UW simulator format
                       % if == 3 print binary
                       % otherwise print both formats
paraview_print = 1;    % if == 1 print paraview visualization
printa         = 1;    % if == 1 save the T matrix = sqrt(lambda)*phi
printabin      = 1;   % if == 1 save the T in a binary file
printa_cond    = 10;   % if == 1 find the best conditioning (PLEASE DO NOT USE)
estatistica    = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(inputbox==1)
    [Lx,Ly,Lz,nx,ny,nz,NX,NY,NZ,eta1,eta2,eta3,...
        varY,Nrand,interpolacao,M,ntipo,TIPOINPUT,...
        file_input_cond]=finputbox();
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(interpolacao~=1)
    NX=nx;
    NY=ny;
    NZ=nz;
    tEINT=0.0;
else
    if(nx==NX)
        nx = nx;
    else
        nx = nx+1;
    end
    if(ny==NY)
        ny = ny;
    else
        ny = ny+1;
    end
    if(nz==NZ)
        nz = nz;
    else
        nz = nz+1;
    end
%
    if(NX<nx)
        disp('Ploblem in the interpolation: NX<nx');
        return
    end
    if(NY<ny)
        disp('Ploblem in the interpolation: NY<ny');
        return
    end
    if(NZ<nz)
        disp('Ploblem in the interpolation: NZ<nz');
        return
    end
end
if(M>nx*ny*nz)
    fprintf('PROBLEM IN M SIZE\n');
    fprintf('ACTUAL M: %d\n',M);
    fprintf('MAXIMUM M: %d\n',nx*ny*nz);
    M = input('Enter new M value: ');
end
if(M<=0)
    M=nx*ny*nz;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define the names of the eigenpairs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ntipo==2)
    coefic=['_b' num2str(beta)];
else
    coefic=['_l' num2str(eta1) 'x' num2str(eta2) 'x' num2str(eta3)];
end
name_aux = [num2str(ntipo,5) '_'...
    num2str(Lx) 'x' num2str(Ly) 'x' num2str(Lz) '_'...
    num2str(NX) 'x' num2str(NY) 'x' num2str(NZ) ...
    coefic '_M' num2str(M) '.dat'];
name_autovet=[homeT 'out/avet' name_aux]
name_autoval=[homeT 'out/aval' name_aux]
name_condP  =[homeT 'out/condP' name_aux]
%**************************************************************************
%** Parameters adjust *****************************************************
sigma    = sqrt(varY);
num_elem = nx*ny*nz;
hx=Lx/double(nx);
hy=Ly/double(ny);
hz=Lz/double(nz);
if(ntipo==2)
    ft = sqrt(2.0)*(varY);
    fat = (varY)*(hx^beta)*((hx/cutoff)^-beta);
else
    fat=1.0;
    ft=varY;
end
%
%**************************************************************************
%**************************************************************************
%** CONDITIONING INPUT ****************************************************
% vet(n,i) coordenada i da posicao do dado condicionado n
if(TIPOINPUT == 1)
    inp = load(file_input_cond);
    vet = inp(:,1:3)
    dados=inp(:,4)
    clear inp
else
    vet=[];
    dados=[];
end
informacoes3D(Lx,Ly,Lz,NX,NY,NZ,vet,dados,beta,eta1,ntipo);
%**************************************************************************
%
%**************************************************************************
%**************************************************************************
%** COVARIANCE MATRIX *****************************************************
%**************************************************************************
%**************************************************************************
disp('------------------------------');
disp('------------------------------');
%disp('CRIANDO A MATRIZ DE COVARIANCIA');
disp('BILDING THE COVARIANCE MATRIX');
tSMatrix=tic;
% constroi a matriz do problema de autovalores 
% Coordinates matrix:
% coord(1:3,j): (x,y,z)- coordinate of the element center
coord=zeros(3,num_elem);
k=0;
for l=1:nz
    for j = 1:ny
        for i = 1:nx
            k = k + 1;
            coord(1,k) = (i-1)*hx+hx/2;
            coord(2,k) = (j-1)*hy+hy/2;
            coord(3,k) = (l-1)*hz+hz/2;
%            fprintf('%d %f %f %f\n',k,coord(1,k),coord(2,k),coord(3,k))
        end
    end
end
wP=1;
%**************************************************************************
% bulding the Covariance matrix *******************************************
C = zeros(num_elem,num_elem);
ephilon = hx;
%
if ntipo==2
    aux1 = (fat*alpha^(beta)*sigma^2*(hx/4)^2*wP');
    for ei = 1:num_elem
        zi = coord(:,ei)';
        C(ei,ei) =  alpha^(beta)*sigma^2*hx^2*sqrt(2);
        xv = [zi(1),zi(2),zi(3)];
        for ej = ei+1:num_elem
            zj = coord(:,ej)';
            yv = [zj(1),zj(2),zj(3)];
            aux = (ephilon./sqrt((xv(:,1)-yv(:,1)).^2+...
                (xv(:,2)-yv(:,2)).^2+(xv(:,3)-yv(:,3)).^2)).^beta;
            C(ei,ej) = aux1*aux;
            C(ej,ei) = C(ei,ej);
        end
    end
end
if ntipo==1
    l1=eta1*eta1;
    l2=eta2*eta2;
    l3=eta3*eta3;
    for ei = 1:num_elem
        zi = coord(:,ei)';
        xv = [zi(1),zi(2),zi(3)];
        for ej = ei:num_elem
            zj = coord(:,ej)';
            yv = [zj(1),zj(2),zj(3)];
             C(ei,ej) = Cov3D(xv,yv,varY,l1,l2,l3);
             C(ej,ei) = C(ei,ej);
        end
    end
end
if ntipo==3
    l1=eta1^2;
    l2=eta2^2;
    l3=eta3^2;
    for ei = 1:num_elem
        zi = coord(:,ei)';
        xv = [zi(1),zi(2),zi(3)];
        for ej = ei:num_elem
            zj = coord(:,ej)';
            yv = [zj(1),zj(2),zj(3)];
            C(ei,ej) = Cov3Df(xv,yv,varY,l1,l2,l3);
            C(ej,ei) = C(ei,ej);
        end
    end
end
if ntipo==4
    l1=eta1*eta1;
    l2=eta2*eta2;
    l3=eta3*eta3;
    aux1=(hx)^2;
    for ei = 1:num_elem
        zi = coord(:,ei)';
        xv = [zi(1),zi(2),zi(3)];
        for ej = ei:num_elem
            zj = coord(:,ej)';
            yv = [zj(1),zj(2),zj(3)];
            C(ei,ej) = aux1*matern(xv,yv,varY,l1,l2,l3,nu);
            C(ej,ei) = C(ei,ej);
        end
    end
end
clear xv yv wP aux coord zi zj
tEMatrix=toc(tSMatrix);
%disp(['Matriz C finalizada em: ' num2str(tEMatrix) ' seg.']);
disp(['C matrix done: ' num2str(tEMatrix) ' seg.']);
disp('------------------------------');
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
disp('------------------------------');
%disp('CALCULO DAS AUTOFUNCOES E AUTOVALORES ');
disp('COMPUTING THE EIGENVALUES AND EIGENVECTORS');
tSauto=tic;
%**************************************************************************
C=(C+conj(C)')/2;
[phi,D] = eig(C);

[lambda,ind] = sort(diag(D),'descend');
phi = phi(:,ind);
clear C D ind

if graf==1
    m = [1:num_elem];
    lambdafig(m,lambda,M);
    if(ntipo==2)
        tipo=[home_fig 'p_autoval'];
        a = num2str(beta,5);
        b=[];
        for ca=1:size(a,2)
            if(a(ca)=='.')
                a;
            else
                b=[b a(ca)];
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(ntipo==1)
        tipo=[home_fig 'e_autoval'];
        a = num2str(eta1,5);
        b=[];
        for ca=1:size(a,2)
            if(a(ca)=='.')
                a;
            else
                b=[b a(ca)];
            end
        end
     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(ntipo==3)
        tipo=[home_fig 'ef_autoval'];
        a = num2str(eta1,5);
        b=[];
        for ca=1:size(a,2)
            if(a(ca)=='.')
                a;
            else
                b=[b a(ca)];
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(ntipo==4)
        tipo=[home_fig 'mf_autoval'];
        a = num2str(eta1,5);
        b=[];
        for ca=1:size(a,2)
            if(a(ca)=='.')
                a;
            else
                b=[b a(ca)];
            end
        end
    end
    name = [tipo num2str(NX,5) 'x' num2str(NY,5) 'x'...
        num2str(NZ,5) '_' b '_' num2str(M,5)];
%    print('-depsc','-r300',name)
end
tEauto=toc(tSauto);
%disp(['Calculo dos auto-valores finalizado em: ' num2str(tEauto) ' seg.']);
disp(['Eigenpairs computation done: ' num2str(tEauto) ' seg.']);
disp('------------------------------');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******** MATRIZ Theta*sqrt(lambda) ***************************************
tMatT=tic;
disp('------------------------------');
%disp('CONSTRUCAO DA MATRIZ T');
disp('BILDING T MATRIX');
T=zeros(num_elem,M);
TOLER = 1.0e-09;
for i=1:M
    if(lambda(i)<TOLER)
        lambda(i) = 0.0;
    end
end
T=phi(:,1:M)*diag(sqrt(lambda(1:M)));
tEMatT=toc(tMatT);
clear phi
%disp([' Matriz T finalizada em: ' num2str(tEMatT) ' seg.']);
disp(['Matrix T done: ' num2str(tEMatT) ' seg.']);
disp('------------------------------');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INTERPOLATION OF T MATRIX TO THE MESH NX x NY x NZ %%%%%%%%%%%%%%%%%%%%
if(interpolacao==1)
    disp('------------------------------');
%    disp('INTERPOLANDO')
    disp('INTERPOLATION')
    tINT=tic;
    if(NX==nx)
        dx = Lx/double(nx);
        x1 = dx/2:dx:Lx;
    else
        dx = Lx/double(nx-1);
        x1 = 0:dx:Lx;
    end
    if(NY==ny)
        dy = Ly/double(ny);
        y1 = dy/2:dy:Ly;
    else
        dy = Ly/double(ny-1);
        y1 = 0:dy:Ly;
    end
    if(NZ==nz)
        dz = Lz/double(nz);
        z1 = dz/2:dz:Lz;
    else
        dz = Lz/double(nz-1);
        z1 = 0:dz:Lz;
    end
%
    [X1,Y1,Z1] = meshgrid(x1,y1,z1);
%
    dx = Lx/double(NX);
    dy = Ly/double(NY);
    dz = Lz/double(NZ);
    x1 = dx/2:dx:Lx;
    y1 = dy/2:dy:Ly;
    z1 = dz/2:dz:Lz;
    [X2,Y2,Z2] = meshgrid(x1,y1,z1);
    clear x1 y1 z1
    newT=zeros(NX*NY*NZ,M);
    for m=1:M
        vect=zeros(ny,nx,nz);
        k=0;
        for l=1:nz
            for j=1:1:ny
                for i=1:nx
                    k=k+1;
                    vect(j,i,l)=T(k,m);
                end
            end
        end
        if(nz==1)
            if(ny==1)
                vect2=interp1(X1,vect,X2,'spline');
                k=0;
                for i=1:NX
                    k=k+1;
                    newT(k,m)=vect2(i);    
                end
            else
                vect2=interp2(X1,Y1,vect,X2,Y2,'spline');
                k=0;
                for j=1:1:NY
                    for i=1:NX
                        k=k+1;
                        newT(k,m)=vect2(j,i);
                    end
                end
            end
        else
            vect2=interp3(X1,Y1,Z1,vect,X2,Y2,Z2,'spline');
            k=0;
            for l=1:NZ
                for j=1:1:NY
                    for i=1:NX
                        k=k+1;
                        newT(k,m)=vect2(j,i,l);
                    end
                end
            end
        end
    end
    tEINT=toc(tINT);
    clear T vect vect2 X1 X2 Y1 Y2 Z1 Z2 x1 y1 z1
    T=newT;
    clear newT
%%% RESCALING THE PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nx=NX;
    ny=NY;
    nz=NZ;
    num_elem = nx*ny*nz;
    disp(['Interpolation done: ' num2str(tEINT) ' seg.']);
    disp('------------------------------');
end
%**************************************************************************
%**************************************************************************
%**** IMPRESSAO DOS AUTOPARES *********************************************
%**************************************************************************
if(printa==1)
%    disp('SALVANDO OS AUTOPARES');
    disp('SAVING THE EIGENPAIRS');
    tSave=tic;
    lambda=lambda(1:M);
    save(name_autoval,'lambda','-ascii');
    if printabin == 1
        name_autovet= [name_autovet(1:end-3) 'bin'];
        fileID=fopen(name_autovet,'w+','l');
        fwrite(fileID,reshape(T,1,num_elem*M),'single');
        fclose(fileID);
    else
        save(name_autovet,'T','-ascii');
    end
    tEsave=toc(tSave);
%    disp(['Autopares salvos em: ' num2str(tEsave) ' seg.']);
    disp(['Eigenpairs saved: ' num2str(tEsave) ' seg.']);
    disp('------------------------------');
else
    tEsave=0.0;
end
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**** LOOP TO REALIZATIONS ************************************************
%**************************************************************************
disp('------------------------------');
%disp('GERACAO DOS CAMPOS');
disp('FILED GENERATION');
tSgera=tic;
n_dados = size(vet,1);
if(ntipo==2)
    if(n_dados>0)
        tipo='pc';
    else
        tipo='p';
    end
%
    a = num2str(beta,5);
    b=[];
    for ca=1:size(a,2)
        if(a(ca)=='.')
            a;
        else
            b=[b a(ca)];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ntipo==1)
    if(n_dados>0)
        tipo='ec';
    else
        tipo='e';
    end
%
    a = num2str(eta1,5);
    b=[];
    for ca=1:size(a,2)
        if(a(ca)=='.')
            a;
        else
            b=[b a(ca)];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ntipo==3)
    if(n_dados>0)
        tipo='efc';
    else
        tipo='ef';
    end
%
    a = num2str(eta1,5);
    b=[];
    for ca=1:size(a,2)
        if(a(ca)=='.')
            a;
        else
            b=[b a(ca)];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ntipo==4)
    if(n_dados>0)
        tipo='mc';
    else
        tipo='m';
    end
%
    a = num2str(eta1,5);
    b=[];
    for ca=1:size(a,2)
        if(a(ca)=='.')
        else
            b=[b a(ca)];
        end
    end
end
g1 = num2str(eta1,5);
g2 = num2str(eta2,5);
g3 = num2str(eta3,5);
if(nz==1)
    name = [tipo num2str(Lx,5) 'x' num2str(Ly,5) '_'...
        num2str(NX,5) 'x' num2str(NY,5)...
        '_l' g1 'x' g2 '_'];
else
    name = [tipo num2str(Lx,5) 'x' num2str(Ly,5) 'x' num2str(Lz,5) '_'...
        num2str(NX,5) 'x' num2str(NY,5) 'x' num2str(NZ,5)...
        '_l' g1 'x' g2 'x' g3 '_'];
end
%
%******** CONDITIONING ****************************************************
% *************************************************************************
% NODAL POSITION OF THE CONDITIONED POINTS ********************************
if(TIPOINPUT == 1)
    n_dados = size(vet,1);
    NT=n_dados;
    A=zeros(n_dados,n_dados);
    vet=vet+1e-6;
    pnode=[];
    for nn=1:n_dados
        k=0;
        for l=1:nz
            zf=l*hz;
            zi=zf-hz;
            if(((vet(nn,3)<zf)&&(vet(nn,3)>zi))||(nz==1))
                for j=ny:-1:1
                    yi = (ny-j)*hy;
                    yf = yi+hy;
                    if((vet(nn,2)<yf)&&(vet(nn,2)>yi))
                        for i=1:nx
                            xf = i*hx;
                            xi = xf - hx;
                            if((vet(nn,1)<xf)&&(vet(nn,1)>xi))
                                k=k+1;
                                pnode(nn)=k;
                            else
                                k=k+1;
                            end
                        end
                    else
                        k=k+nx;
                    end
                end
            else
                k=k+nx*ny;
            end
        end
    end
    pnode = pnode';
%%
% *** vetor de posiccoes em relaccao aos nos da malha *********************
% *** apenas para nodes ***************************************************
    pt = zeros(num_elem-n_dados,1);
    k=0;
    for i=1:num_elem
        sgn = 1;
        for j=1:n_dados
            if (i == pnode(j))
                sgn = -1;
                break;
            end
        end
        if sgn>0
            k=k+1;
            pt(k)=i;
        end
    end
end
%**************************************************************************
%******** LOOP ************************************************************
    mu=0.0;
    sig=1.0;
    Xi=zeros(num_elem,1);
    auxY = T(1,1:M)*T(1,1:M)'/2 + sigma^2/2;
    mY=zeros(num_elem,1);
    X = [];
%**************************************************************************
%** PROCURAR O MELHOR CONDICIONAMENTO *************************************
if(n_dados>0)
    if(printa_cond==1)
        [vect,vect2,tCOND]=bestMatCond(M,NT,T,pnode);
        fid = fopen(name_condP, 'wt');
        for i=1:NT
            fprintf(fid, '%d\n', vect(i));
        end
        for i=1:M-NT
            fprintf(fid, '%d\n', vect2(i));
        end
        fclose(fid);
    else
        tCoMatrix=tic;
        vect =[1:n_dados];
        vect2=[1+n_dados:M];
        tCOND=toc(tCoMatrix);
    end
else
    tCOND=0.0;
end
%**************************************************************************
%**************************************************************************
for nr=1:Nrand
%******** CONDICIONAMENTO DO CAMPO ****************************************
        theta = lhsnorm(mu,sig,M);
        if(n_dados>0)
            b=zeros(n_dados,1);
            for i=1:n_dados
                for j=1:M-NT
                  b(i)= b(i)-T(pnode(i),vect2(j))*theta(vect2(j));
               end
               b(i)=b(i)+dados(i)+auxY-sigma^2/2;
               for j=1:NT
                  A(i,j)=T(pnode(i),vect(j));
               end
            end
            alpha=inv(A)*b;
            for i=1:NT
               theta(vect(i))=alpha(i);
            end
        end
%**************************************************************************
        for el=1:num_elem
            Xi(el) = mY(el) + T(el,:)*theta;
            %Xi(el) = mY(el) + T(el,:)*theta - auxY;
        end
        if(estatistica==1)
            MEDIA = mean(Xi)
            VAR   = var(Xi)
            X = [X; Xi];
        end
%******* PRINTING THE FIELDS **********************************************
        if(nz==1)
            imprime(Lx,Ly,nx,ny,ntipo,beta,Xi,nr,home,name,tipo_prt);
        else
            imprime3D(Lx,Ly,Lz,nx,ny,nz,ntipo,beta,Xi,nr,home,name,tipo_prt);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAVIEW PRINTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(paraview_print==1)
            filename = [homep 'field_' name 'M' num2str(M,5)...
                '-' num2str(nr-1,5)];
            paraviewprinter(Lx,Ly,Lz,nx,ny,nz,Xi,filename)
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
clear T;
tEgera=toc(tSgera);
tElapsed=toc(tStart);
disp(['End of fields generation: ' num2str(tEgera) ' seg.'])
disp('------------------------------');

if(estatistica==1)
    fprintf('Mean of Y.......: %f\n',mean(X));
    fprintf('Variance of Y...: %f\n',var(X));
    fprintf('Minimum of Y....: %f\n',min(X));
    fprintf('Maximum of Y....: %f\n',max(X));
end

if(fig==1&&nz~=1)
    dx = Lx/double(NX);
    dy = Ly/double(NY);
    dz = Lz/double(NZ);
    x1 = dx/2:dx:Lx;
    y1 = dy/2:dy:Ly;
    z1 = dz/2:dz:Lz;
    [X2,Y2,Z2] = meshgrid(x1,y1,z1);
    X=X2*0;
    k=0;
    for l=1:NZ
        for j = 1:NY
            for i = 1:NX
                k=k+1;
                X(j,i,l)=Xi(k);
            end
        end
    end
    yslice = [Ly/2];
    zslice = [Lz/2];
    xslice = [dx, Lx/2, Lx-dx];
    figure(2)
    slice(X2,Y2,Z2,X,xslice,yslice,zslice), shading flat;
    daspect([ 1 1 1]);
    xlim([0 Lx]);
    ylim([0 Ly]);
    zlim([0 Lz]);
    view(-35,20);
    name = [home_fig tipo '_field_' num2str(NX,5),...
        'x' num2str(NY,5) 'x' num2str(NZ,5) '_' num2str(M,5)];
    print('-depsc','-r300',name);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('############################################');
%disp(['TEMPO TOTAL GASTO: ' num2str(tElapsed) ' seg.']);
disp(['TOTAL TIME: ' num2str(tElapsed) ' seg.']);
disp('############################################');
%disp(['PORCENTAGEM DO TEMPO TOTAL GASTO PARA CONSTRUIR A MATRIZ DE COVARIANCIA: '...
%    num2str(100*tEMatrix/tElapsed) ' %']);
disp(['PERCENTAGE OF TOTAL TIME SPENT TO CONSTRUCTION OF THE C MATRIX: '...
    num2str(100*tEMatrix/tElapsed) ' %']);
disp('############################################');
%disp(['PORCENTAGEM DO TEMPO TOTAL GASTO PARA O CALCULO DOS AUTOPARES: '...
%    num2str(100*tEauto/tElapsed) ' %']);
disp(['PERCENTAGE OF TOTAL TIME SPENT TO COMPUTE THE EIGENPAIRS: '...
    num2str(100*tEauto/tElapsed) ' %']);
disp('############################################');
%disp(['PORCENTAGEM DO TEMPO TOTAL GASTO PARA CONSTRUIR A MATRIZ T (phi*sqrt(lambda)): '...
%    num2str(100*tEMatT/tElapsed) ' %']);
disp(['PERCENTAGE OF TOTAL TIME SPENT TO CONSTRUCTION OF THE T MATRIZ (phi*sqrt(lambda)): '...
    num2str(100*tEMatT/tElapsed) ' %']);
disp('############################################');
%disp(['PORCENTAGEM DO TEMPO TOTAL GASTO PARA INTERPOLACAO DA MATRIZ T PARA MALHA FINA: '...
%    num2str(100*tEINT/tElapsed) ' %']);
disp(['PERCENTAGE OF TOTAL TIME SPENT TO INTERPOLATION: '...
    num2str(100*tEINT/tElapsed) ' %']);
disp('############################################');
%disp(['PORCENTAGEM DO TEMPO TOTAL GASTO PARA ESCOLHER O MELHOR CONDICIONAMENTO MATRIZ: '...
%    num2str(100*tCOND/tElapsed) ' %']);
%disp('############################################');
%disp(['PORCENTAGEM DO TEMPO TOTAL GASTO PARA SALVAR A MATRIZ T: '...
%    num2str(100*tEsave/tElapsed) ' %']);
disp(['PERCENTAGE OF TOTAL TIME SPENT TO SAVE THE T MATRIX: '...
    num2str(100*tEsave/tElapsed) ' %']);
disp('############################################');
%disp(['PORCENTAGEM DO TEMPO TOTAL GASTO PARA GERAR OS CAMPOS: '...
%    num2str(100*tEgera/tElapsed) ' %']);
disp(['PERCENTAGE OF THE TOTAL TIME SPENT TO GENERATE THE FIELDS: '...
    num2str(100*tEgera/tElapsed) ' %']);
disp('############################################');
%**************************************************************************
clear;
