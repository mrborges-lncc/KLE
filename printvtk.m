function printvtk(nos,pos,nome,X)
  file = nome;
  fid = fopen(file,'w');
  fprintf(fid,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n');
  fprintf(fid,'POINTS %10d float\n',size(pos,1));
  for i=1:size(pos,1)
      fprintf(fid,' %10.8e  %10.8e  %10.8e\n',pos(i,1),pos(i,2),pos(i,3));
  end
  fprintf(fid,'CELLS %10d %10d\n', size(nos,1),size(nos,1)*(8+1));
  for i=1:size(nos,1)
      fprintf(fid,' %10d %10d %10d %10d %10d %10d %10d %10d %10d\n',8,nos(i,1)-1,nos(i,2)-1,...
          nos(i,3)-1,nos(i,4)-1,nos(i,5)-1,nos(i,6)-1,nos(i,7)-1,nos(i,8)-1);
  end
  fprintf(fid,'CELL_TYPES %10d \n', size(nos,1));
  for i=1:size(nos,1)
      fprintf(fid,' %d\n',12);
  end
  fprintf(fid,'CELL_DATA %10d\n', size(nos,1));
%  fprintf(fid,'SCALARS log(K) float\n');
  fprintf(fid,'SCALARS log(K) double\n');
  fprintf(fid,'LOOKUP_TABLE default\n');

    for i=1:size(nos,1)
        fprintf(fid,'%g\n',X(i));
     end

%   for i=1:size(nos,1)
%       fprintf(fid,' %f\n',i);
%   end
  fclose(fid);