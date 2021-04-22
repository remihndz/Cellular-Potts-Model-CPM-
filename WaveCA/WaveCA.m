% Creates an animation of wave propagation using a simple cellular
% automaton (CA). 
% The algorithm is described in: 
% Komatsuzaki T., Iwata Y., Morishita S. (2012) Modelling of Incident Sound
% Wave Propagation around Sound Barriers Using Cellular Automata.In: 
% Sirakoulis G.C., Bandini S. (eds) Cellular Automata. ACRI 2012. Lecture
% Notes in Computer Science, vol 7495. Springer, Berlin, Heidelberg. 
% https://doi.org/10.1007/978-3-642-33350-7_40
clear all;

%% Parameters and initialisation
damping = 0.001;
ca      = 1/sqrt(2);
omega   = 3/(2*pi);
dt      = 0.3;

Max_Steps = 1000;
freq_plot = 10;

% Number of cells in the extended lattice
Width = 10; Height = 8; % in meters
Scale = 100; % i.e. Scale #cells = 1m
Extra_cells = 2*Scale; % Number of extra cells to simulate inifnite size domain

Rows = Scale*Height; 
Cols = Scale*Width;

% Size of the analysed part of the lattice
Row_0 = Extra_cells; Row_1 = Rows-Extra_cells;  
Col_0 = Extra_cells; Col_1 = Cols-Extra_cells;

Pressure = zeros(Rows,Cols);
Velocity = zeros(Rows, Cols, 4);
[Source_x, Source_y] = deal(Row_0, Col_0);

im = imagesc(Pressure(Row_0:Row_1, Col_0:Col_1));
colorbar
caxis([0,0.2])
colormap('jet')
drawnow
hold on


%% Position of accoustic walls
% Main branch
Wall_x = Source_x + 100; Wall_y = Source_y + 100;
Wall_length = 2*Scale;
Wall = ones(Rows,Cols);
Wall(Wall_y:Wall_y+Wall_length, Wall_x) = 0.0;

Branch_length = floor(Wall_length/2);
for i=1:Branch_length
    Wall(Wall_y-i, Wall_x+i) = 0.0;
    Wall(Wall_y-i, Wall_x-i) = 0.0;
end

SubBranch_length = floor(Branch_length/4);

Start_SubBranch1_x = Wall_x-2*floor(Branch_length/3);
Start_SubBranch1_y = Wall_y-2*floor(Branch_length/3);

Start_SubBranch2_x = Wall_x+2*floor(Branch_length/3);
Start_SubBranch2_y = Wall_y-2*floor(Branch_length/3);

for i=1:SubBranch_length
    Wall(Start_SubBranch1_y-i, Start_SubBranch1_x+i) = 0.0; 
    Wall(Start_SubBranch2_y-i, Start_SubBranch2_x-i) = 0.0; 
end


% End points of wall segments (for plotting)
Segment_x = [];
Segment_y = [];

Segment_x = [Segment_x Wall_x Wall_x];
Segment_y = [Segment_y Wall_y Wall_y+Wall_length];

Segment_x = [Segment_x Wall_x Wall_x+Branch_length];
Segment_y = [Segment_y Wall_y Wall_y-Branch_length];

Segment_x = [Segment_x Wall_x Wall_x-Branch_length];
Segment_y = [Segment_y Wall_y Wall_y-Branch_length];

Segment_x = [Segment_x Start_SubBranch1_x Start_SubBranch1_x+SubBranch_length];
Segment_y = [Segment_y Start_SubBranch1_y Start_SubBranch1_y-SubBranch_length];

Segment_x = [Segment_x Start_SubBranch2_x Start_SubBranch2_x-SubBranch_length];
Segment_y = [Segment_y Start_SubBranch2_y Start_SubBranch2_y-SubBranch_length];

Segment_x = Segment_x-Extra_cells; Segment_y = Segment_y-Extra_cells;

%% Simulation
for k=1:Max_Steps

    Pressure(Source_x, Source_y) = sin(omega*k*dt);
    
    % Update outflow velocities
    % Upward flow
    Velocity(:, 1:Cols-1, 1) = Velocity(:, 1:Cols-1, 1) ...
        + Pressure(:,1:Cols-1) - Pressure(:,2:Cols);
    Velocity(:, Cols, 1) = Pressure(:, Cols);
    
    % Downward flow
    Velocity(:, 2:Cols, 2)  = Velocity(:, 2:Cols, 2)  ...
        + Pressure(:,2:Cols) - Pressure(:,1:Cols-1);
    Velocity(:, 1, 2) = Pressure(:, 1);
    
    % Leftward flow
    Velocity(2:Rows, :, 3)  = Velocity(2:Rows, :, 3)  ...
        + Pressure(2:Rows, :) - Pressure(1:Rows-1,:);
    Velocity(1, :, 3) = Pressure(1, :);
    
    
    % Rightward flow
    Velocity(1:Rows-1, :, 4) = Velocity(1:Rows-1, :, 4) ...
        + Pressure(1:Rows-1, :) - Pressure(2:Rows, :);
    Velocity(Rows, :, 4) = Pressure(Rows, :);
    
    for i = 1:4
        Velocity(:,:,i) = Velocity(:,:,i).*Wall;
    end
    
    % Linear energy dissipation
    Velocity = (1-damping).*Velocity;
    
    % Update Pressure according to directional velocities.
     Pressure = Pressure - ca^2*sum(Velocity,3);
    
    % Update plot
    if rem(k, freq_plot)==0
        im.CData = Pressure(Row_0:Row_1, Col_0:Col_1);
        line(Segment_x(1:6), Segment_y(1:6), 'Color', 'k', 'LineWidth', 1.5);
        line(Segment_x(7:8), Segment_y(7:8), 'Color', 'k', 'LineWidth', 1.5);
        line(Segment_x(9:10), Segment_y(9:10), 'Color', 'k', 'LineWidth', 1.5);
        drawnow
    end
end
