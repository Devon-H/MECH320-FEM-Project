%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mario Bras
% March 2023
%
% This code is to be used in the MECH 320 course at UVic
% and may not be shared, uploaded, or distributed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear command window and variables, and close all open figures
clc
clear all
close all

%% 1: Model parameters
half_crack_height = 6.5;       % Half crack height (mm)
half_crack_width = 1;        % Half crack width (mm)
half_width = 38/2;              % Plate half width (mm)
half_length = 240/2;             % Plate half length (mm)
youngs_modulus = 200*10^3;          % Youngs Modulus of material (MPa)
poissons_ratio = 0.3;          % Poissons ration of material

surf_traction = 10;           % Surface traction (N/mm)
mesh_size = half_crack_height/6; % Mesh refinement (no need to change)

plot_scale_factor = 20;       % Scale factor for plotting deformations, adjust it to be able to visualize contraction of the plate near the hole

%% 2: Geometry definition
% See https://www.mathworks.com/help/pde/ug/create-geometry-at-the-command-line.html for more information
% on how to use a geometry description matrix (GDM) to represent the geometry of the problem

% Vector representing rectangular plate geometry
plate = [3 4 -half_length half_length half_length -half_length -half_width -half_width half_width half_width]';

% Vector representing circular or elliptical hole geometry
hole = [1 0 0 half_crack_height]';
hole = [hole;zeros(length(plate)-length(hole),1)];

% Define geometry description matrix (GDM) as a subtraction operation between hole and bar
gdm = [plate, hole];                   % Form the GDM containing all the shapes
names = char('bar','hole')';          % Give names to shapes
geom = decsg(gdm,'bar - hole',names); % Boolean subtraction of hole from bar

%% 3: FEM model (no need to modify this section)
model = createpde('structural','static-planestress');

% Create geometry and include it into the structural model
geometryFromEdges(model,geom);

% Specify material properties
structuralProperties(model,'YoungsModulus',youngs_modulus,'PoissonsRatio',poissons_ratio);

% Apply boundary conditions
structuralBC(model,'Edge',3,'XDisplacement',0);
structuralBC(model,'Vertex',3,'YDisplacement',0);
structuralBoundaryLoad(model,'Edge',1,'SurfaceTraction',[surf_traction;0]);

% Generate mesh
generateMesh(model,'Hmax',mesh_size);

% Solve system of equations
sol = solve(model);

%% 4: Plotting

% Notes:
% - sol.Stress is a structure that contains stress information for the model (e.g. you can access
%   normal and shear stresses from the sol.Stress.sxx, sol.Stress.syy and sol.Stress.sxy fields.
% - sol.Displacement is a structure that contains displacement information for this model. Inspect
%   this structure in the MATLAB workspace to see what fields are available inside it.

% Geometry
figure           % Create figure
pdegplot(model); % Plot model geometry
axis equal       % Set equal scales for x and y axes
title('Plate Geometry');    % Plot title

% Mesh
figure          % Create figure
pdemesh(model); % Plot model mesh
axis equal      % Set equal scales for x and y axes
title('Plate Mesh');   % Plot title

% Von Mises stress contours
svon = 0.5.*sqrt((sol.Stress.sxx-sol.Stress.syy).^2+(0-sol.Stress.sxx).^2+6*(sol.Stress.sxy).^2);   % Calculate von Mises stress

figure                                         % Create figure
pdeplot(model,'XYData',svon,'ColorMap','jet'); % Plot xy data (in this case the von Mises stress)
axis equal                                     % Set equal scales for x and y axes
title('Von Mises Stress');                     % Plot title

% Deformed shape
figure                                         % Create figure
pdeplot(model,"XYData",sol.Displacement.x,"Deformation",sol.Displacement,"DeformationScaleFactor",plot_scale_factor,"ColorMap","jet"); % Plot xy data (in this case displacement and deformation superimposed)
axis equal                                     % Set equal scales for x and y axes
title('Plate Deformation');                    % Plot title

%% 5: Calculate stress concentrantion factor
p = [-100,0];                         % Coordinates of point far away from the hole
snom = interpolateStress(sol,p').xx; % Nominal stress x obtained at point far away from hole
smax = max(sol.Stress.sxx);        % Maximum normal stress x in the plate
%k = 1+2*(half_crack_height/half_crack_width);     % Stress concentration factor for ellipse
k = 1+2*(1);
fprintf("k = %d\n",k);               % Print k to the console
