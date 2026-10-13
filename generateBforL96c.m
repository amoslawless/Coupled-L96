function [Bc] = generateBforL96c(Nx,Ls,variance)
%   Construct climatological Bc from exponential kernal functionmax
%   Corrected by ASL for periodic domain 25/10/24
    xmin = 0; xmax = Nx-1; dx = 1;
    x_grid = xmin:dx:xmax+1;             
    Bc = zeros(Nx, Nx);
    a = Nx*dx;
    dtheta = (2*pi)/Nx;

    % for j =1:Nx
    %     for i = 1:Nx
    %         rad = (abs(x_grid(i)-x_grid(j)));
    %         rad = min(rad,a-rad);        % Periodic domain
    %         Bc(i,j) = variance * exp(-rad/Ls);
    %     end
    % end
    % 
    for i=1:Nx
        for j=1:Nx
            theta(i,j) = abs(i-j) * dtheta;
        end
    end

    for i=1:Nx
        for j=1:Nx
            Bc(i,j) = variance * ( 1 +  abs(2*a*sin(theta(i,j)/2))  / Ls ) ...
                * exp (   -  abs(2*a*sin(theta(i,j)/2))  / Ls );
        end
    end

end

