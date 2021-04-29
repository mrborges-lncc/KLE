function paraviewprinter(Lx,Ly,Lz,nx,ny,nz,X,arq)
%
    numel=nx*ny*nz;
    numno=(nx+1)*(ny+1)*(nz+1);
    bb2 = [arq '.vtk']
%     outfile2 = fopen(bb2, 'wt');
%
%     fprintf(outfile2,'# vtk DataFile Version 3.0\n');
%     fprintf(outfile2,'vtk output\n');
%     fprintf(outfile2,'ASCII\n');
%     fprintf(outfile2,'DATASET UNSTRUCTURED_GRID\n');
%     fprintf(outfile2,'POINTS %d float\n',numno);
% Matriz de coordenadas:
% coord(): (x,y,z)- coordenadas dos centros dos elementos
disp('-------------------------------------------------------------------')
tic
hx=Lx/nx;
hy=Ly/ny;
hz=Lz/nz;
hx2=hx*0.5;
hy2=hy*0.5;
hz2=hz*0.5;
coord=zeros(numno,3);
k = 0;
for z = 1:nz+1
  for j = 1:ny+1
    for i = 1:nx+1
      k = k + 1;
      coord(k,1) = (i-1)*hx;
      coord(k,2) = (j-1)*hy;
      coord(k,3) = (z-1)*hz;
    end
  end
end
%      
%     imprime as coordenadas
%
% l=0;
% for k=1:nz+1
%     for j=1:ny+1
%         for i=1:nx+1
%             %l=i+nx*(j-1)+nx*ny*(k-1);
%             l=l+1;
%             xx=coord(1,l);
%             yy=coord(2,l);
%             zz=coord(3,l);
%             fprintf(outfile2,'%12.5e%12.5e%12.5e\n',xx,yy,zz);
%         end
%     end
% end
% %
nen = 8;
% fprintf(outfile2,'CELLS %8d %8d\n',numel,numel*(nen+1));
% %       
nos = zeros(numel,nen);
l=0;
for k=1:nz
    for j=1:ny
        for i=1:nx
            %l =i+nx*(j-1)+nx*ny*(k-1);
            l=l+1;
            n1=i+(j-1)*(nx+1)+(k-1)*(nx+1)*(ny+1);
            n2=n1+1;
            n3=n2+(nx+1);
            n4=n3-1;
            n5=n1+(nx+1)*(ny+1);
            n6=n5+1;
            n7=n6+(nx+1);
            n8=n7-1;
            nos(l,:)=[n1 n2 n3 n4 n5 n6 n7 n8];
%             fprintf(outfile2,'%8d%8d%8d%8d%8d%8d%8d%8d%8d\n',...
%                 nen,n1-1,n2-1,n3-1,n4-1,n5-1,n6-1,n7-1,n8-1);
        end
    end
end
%
% fprintf(outfile2,'CELLS_TYPES %8d\n',numel);
% for k=1:numel
%     fprintf(outfile2,'%d\n',12)
% end
% %
% fprintf(outfile2,'CELLS_DATA %8d\n',numel);
% fprintf(outfile2,'SCALARS Y float\n');
% fprintf(outfile2,'LOOKUP_TABLE default\n');
% for k=1:numel
%     fprintf(outfile2,'%f\n',k);
% end
% fclose(outfile2);
printvtk(nos,coord,bb2,X)