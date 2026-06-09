clc
clear
close all
format long

%Gitrm History file
path1="/lore/nathd2/DIII-D-helicon/Output/";
file1="gitrm-history.nc";
fid=path1+file1;
x = ncread(fid, 'x');
y = ncread(fid, 'y');
z = ncread(fid, 'z');
sp=ncread(fid, 'species');

nP=size(x,2);
nT=size(x,1);
x1=reshape(x,[1, nP*nT]);
y1=reshape(y,[1, nP*nT]);
z1=reshape(z,[1, nP*nT]);
sp1=reshape(sp,[1, nP*nT]);

%% Gitrm Initial file to know which particle has what mass
mass_id=sp1(1,1:nP*(nT-1));


%% Create Vtk file for visualization
n_cells=nP*(nT-1);
cells=zeros(n_cells,3);
cell_type=zeros(n_cells,1);

cumul=0;
for i=1:nP
  for j=1:nT
    if j<=nT-1
      index=(i-1)*(nT-1)+j;
      cells(index,:)=[2 cumul cumul+1];
      cell_type(index)=3 ;
    end
    cumul=cumul+1 ;
  end
end
fileID = fopen(path1 + "XX1.vtk",'w');

fprintf(fileID,'%6s\n','# vtk DataFile Version 2.0');
fprintf(fileID,'%6s\n','particlepaths');
fprintf(fileID,'%4s\n','ASCII');
fprintf(fileID,'%4s\n','DATASET UNSTRUCTURED_GRID');

fprintf(fileID,'%4s %d %4s \n','POINTS', nP*nT, 'double');
fprintf(fileID,'%0.15f %0.15f %0.15f\n',[x1;y1;z1]);

fprintf(fileID,'%4s %d %d \n','CELLS', n_cells, 3*n_cells);
fprintf(fileID,'%d %d %d \n', cells');
fprintf(fileID,'%4s %d\n','CELL_TYPES', n_cells);
fprintf(fileID,'%d\n',cell_type);

fprintf(fileID,'%4s %d \n','CELL_DATA', nP*(nT-1));
fprintf(fileID,'SCALARS mass int 1\n');
fprintf(fileID,'LOOKUP_TABLE default\n');
fprintf(fileID,'%d\n',mass_id);

fclose(fileID);