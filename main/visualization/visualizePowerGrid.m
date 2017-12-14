function visualizePowerGrid(prefix, mpc, busLocations, LineVars, LineVarsLabels, BusVars, BusVarsLabels)

PointList        = busLocations;
ElementList      = mpc.branch(:, 1:2);
numberOfPoints   = size(busLocations, 1);
numberOfElements = size(ElementList, 1);
numberOfVertices = 2;
type             = 3;  % line segments

Indices          = busLocations(:,1);
nb               = max(Indices);
I                = zeros(nb,1);
I(Indices)       = 1:length(Indices);


%  write the VTK header
% -----------------------
filename = strcat(prefix, '.vtk');
fprintf('writting mesh file %s\n', filename);

fout = fopen(filename, 'w');
fprintf(fout,'# vtk DataFile Version 5.10\n');
fprintf(fout,'Hexahedral mesh with data\n');
fprintf(fout,'ASCII\n');
fprintf(fout,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fout,'POINTS %d float\n', numberOfPoints);


%  write the PointList (x,y coordinates of the buses)
% ---------------------------
for i = 1:numberOfPoints
    x_i = PointList(i,2:3); %bus [x, y]
    fprintf(fout,'%25.16e %25.16e %25.16e\n', x_i(1), x_i(2), 0.0);
end

%  write the Elements (a.k.a branches: from, to)
% -----------------------------------------------
fprintf(fout,'\n');
entries = (numberOfVertices+1)*numberOfElements;
fprintf(fout,'CELLS %d %d\n', numberOfElements, entries);
first_number = 1;

for e = 1:numberOfElements
    
    v_e = ElementList(e, :);
    
    v_e = I(v_e) - first_number;
    
    fprintf(fout,'%d ', numberOfVertices);
    
    for i=1:numberOfVertices
        fprintf(fout,'%d ', v_e(i));
    end
    
    fprintf(fout, '\n');
    
end

%  write the type of elements (type=3 line segment)
% --------------------------------------------------
fprintf(fout,'\n');
fprintf(fout,'CELL_TYPES %d\n', numberOfElements);

for e = 1:numberOfElements
    
    fprintf(fout,'%d\n', type);
    
end
fprintf(fout,'\n');


%  write line properties
% -----------------------
nfields = size(LineVars, 2);

fprintf(fout,'CELL_DATA %d\n', size(LineVars,1)); %numberOfElements);
for f = 1:nfields
  fprintf(fout,'SCALARS %s float 1\n', LineVarsLabels{f});
  fprintf(fout,'LOOKUP_TABLE default\n');
  fprintf( fout,'%25.16e\n', LineVars(:, f) );
end


%  write node properties
% -----------------------
nfields = size(BusVars, 2);
fprintf(fout,'POINT_DATA %d\n', size(BusVars,1)); %numberOfPoints);
for f = 1:nfields
  fprintf(fout,'SCALARS %s float 1\n', BusVarsLabels{f});
  fprintf(fout,'LOOKUP_TABLE default\n');
  fprintf( fout,'%25.16e\n', BusVars(:, f) );
end

fclose(fout);
end

