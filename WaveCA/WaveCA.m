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
damping = 0.9;
omega   = 3/(2*pi);
dt      = 0.3;

Max_Steps = 500;
freq_plot = 5;
[Rows,Cols] = deal(400,400); % Number of cells in the lattice

Pressure = zeros(Rows,Cols); % Average atmospheric pressure
Pressure(:,[1,Cols]) = 0.0;
Pressure([1,Rows],:) = 0.0;
Velocity = zeros(Rows, Cols, 4);
[Source_x, Source_y] = deal(floor(Rows/3), floor(Cols/2));

im = imagesc(Pressure);
colorbar
caxis([-1,1])
colormap('bone')
drawnow

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
    
    % Update Pressure according to directional velocities.
    Pressure = Pressure - 0.5*damping*sum(Velocity,3);
    
    % Update plot
    if rem(k, freq_plot)==0
        im.CData = Pressure;
        drawnow
    end
end
