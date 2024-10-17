function z = S3(x)
    % Define constants
    g = 9.81;

    % Define initial conditions for differential equations.
    x0 = x;
    y0 = 2 * x0 .* exp(-8 * x0.^4) + x0.^2 / 2;
    % Set initial velocity to zero
    vx0 = 0;
    vy0 = 0;

    % Create initial condition vector.
    z0 = [x0, y0, vx0, vy0];

    % Solve equations from t = 0 to t = 100.
    t = linspace(0, 1, 101);

    % Set accuray of ode45.
    opts = odeset('RelTol', 1e-12);

    % Solve differential equation.
    [t, z] = ode45(@rhs, t, z0, opts);

    function [dz] = rhs(~, z)
        % Seperate z
        x = z(1);
        vx = z(3);
        vy = z(4);
        
        dy = exp(-8*x^4) * (x*exp(8*    x^4) - 64*x^4 + 2);
        d2y = exp(-8*x^4) * (exp(8*x^4) + 2048*x^7 - 320*x^3);
        
        % Compute accelerations in x and y 
        ax = -((d2y .* vx.^2 + g) .* dy) ./ (1 + dy.^2);
        ay =  d2y .* vx.^2 + dy .* ax;
        
        % Return derivatives
        dz = zeros(4,1);
        dz(1) = vx;
        dz(2) = vy;
        dz(3) = ax;
        dz(4) = ay;
    end

    % Create figure for animation.
    figure;
    hold on;
    axis equal;
    xlim([-5, 5]);
    ylim([-1, 15]);
    grid on;
    xlabel('X Position')
    ylabel('Y Position')

    % Plot wire shape
    x_vec = linspace(-10,10,10000);
    y_vec = 2 * x_vec .* exp(-8 * x_vec.^4) + x_vec.^2 / 2;
    plot(x_vec,y_vec)
    
    % Initilize animations.
    bead = plot(z(1,1), z(1,2), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    g_vec_plot = quiver(z(1,1), z(1,2), 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 2);
    normal_force_plot = quiver(z(1,1), z(1,2), 0, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 2);
    net_force_plot = quiver(z(1,1), z(1,2), 0, 0, 'm', 'LineWidth', 2, 'MaxHeadSize', 2);

    legend_entries = [g_vec_plot, normal_force_plot, net_force_plot];
    legend_labels = {'Gravitational Force', 'Normal Force', 'Net Force'};
    legend(legend_entries, legend_labels, 'Location', 'southoutside');
    
    % Animate
    for i = 1:length(t)
        set(bead, 'XData', z(i,1), 'YData', z(i,2))
        set(g_vec_plot, 'XData', z(i,1), 'YData', z(i,2), 'UData', 0, 'VData', -g/g);
        
        % Calculate normal vector components
        dz = rhs(t(i),z(i,:));

        Nx = dz(3); 
        Ny = dz(4) + g;      

        % Normalize the normal vector
        Nx = Nx / g;
        Ny = Ny / g;

        % Update normal force vector
        set(normal_force_plot, 'XData', z(i,1), 'YData', z(i,2), 'UData', Nx, 'VData', Ny);
        
        % Calculate net force 
        F_net_x = Nx;           
        F_net_y = Ny - g/g;      
        set(net_force_plot, 'XData', z(i,1), 'YData', z(i,2), 'UData', F_net_x, 'VData', F_net_y);
        
        title_name = sprintf('Bead on Wire Simulation: Time = %.2f s', t(i));
        title(title_name);
        drawnow;
    end
end