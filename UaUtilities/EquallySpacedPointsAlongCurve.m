

function [pt_x,pt_y,normal_x,normal_y]=EquallySpacedPointsAlongCurve(x,y,num_points)


x=x(:);
y=y(:);


% 3. Calculate cumulative arc length 
distances = [0; cumsum(sqrt(diff(x).^2 + diff(y).^2))];

% 4. Interpolate evenly spaced distances
target_distances = linspace(0, distances(end), num_points);
pt_x = interp1(distances, x, target_distances, 'linear');
pt_y = interp1(distances, y, target_distances, 'linear');

% 5. Calculate the unit tangent (dx, dy) and rotate to get the normal (-dy, dx)
dx = gradient(pt_x);
dy = gradient(pt_y);
tangent_magnitudes = sqrt(dx.^2 + dy.^2);
unit_tangent_x = dx ./ tangent_magnitudes;
unit_tangent_y = dy ./ tangent_magnitudes;

% Normal vectors (rotated 90 degrees counter-clockwise)
normal_x = -unit_tangent_y;
normal_y = unit_tangent_x;


return

% 6. Plotting
figure(10);
plot(x, y, 'b', 'LineWidth', 1.5); hold on;
scatter(pt_x, pt_y, 40, 'r', 'filled'); % Equidistant points

% Scale the normal vectors for visualization
scale = 0.5; 
quiver(pt_x, pt_y, normal_x * scale, normal_y * scale, 0, 'Color', 'g', 'LineWidth', 1.5);

title('Equally Distributed Points and Normals');
xlabel('X'); ylabel('Y');
legend('Original Curve', 'Equidistant Points', 'Normals');
axis equal; grid on;
hold off;

end