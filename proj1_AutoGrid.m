function [condstested,result ] = proj1_AutoGrid(y3_0,y4_0, alpha, beta, nu,kappa, NumSteps,LenGrid, scale, gridorigin_x, gridorigin_y)
% Runs proj1_func.m for many sets of initial conditions, and records
% whether they are stable or unstable.
% We fix alpha, beta, nu, and y3_0, y4_0, and this checks a grid of varying
% y1_0 and y2_0 initial conditions, testing a 2d grid in 4d space
% If one wishes to fix a different set of two variables to vary, one can
% manually change the function.
% To determine what grid of conditions to check, gridorigin_x and _y
% determine the center of this grid, LenGrid determines the number of
% points in each side of the square lattice, and "scale" is a positive
% number that linearly scales the size of the grid (i.e. a scale of 1 sets
% the lattice to have a distance of 1 between points, and a scale of .5
% sets a distance of .5 between each point)
% Example 
% proj1_AutoGrid(0, 0, 1, 1, 1, 0, 500, 4,0.5,3,0);


initconds = zeros(4,LenGrid*LenGrid);
status = true(1,length(initconds));

%Fills in the grid of initial conditions as specified above
for i = 1:LenGrid 
    for j = 1:LenGrid
        initconds(1,(i-1)*LenGrid+j) = (i-(LenGrid+1)/2)*scale+gridorigin_x;
        initconds(2,(i-1)*LenGrid+j) = (j-(LenGrid+1)/2)*scale+gridorigin_y;
        initconds(3,(i-1)*LenGrid+j) = y3_0;
        initconds(4,(i-1)*LenGrid+j) = y4_0;
    end
end


for k = 1:length(initconds) %Calculates maximum radius each trajectory acheives
    crdn = proj1_func(initconds(1,k),initconds(2,k),initconds(3,k),initconds(4,k),alpha,beta,nu,kappa,NumSteps);
    M = sqrt(crdn(1,length(crdn))^2+crdn(2,length(crdn))^2);
    
    % Check if each final distance from origin is past the limit (which means
    % it explodes)
    if M < 10^4 %picked radius of 200 as Exloding criterion
        status(k) = true;
    else
        status(k) = false;
    end 
    
end

if all(status) %Determines if any trajectories exploded
    result = 'Good';
    figure(2)
    hold on
    plot(initconds(1,:),initconds(2,:),'b*')
    xlabel('$y_1$','FontSize', 12, 'interpreter','latex') % x-axis label
    ylabel('$y_2$','FontSize', 12, 'interpreter','latex') % y-axis label
    legend('Stable');
    hold off
    condstested = [initconds(1,:);initconds(2,:)]';
else 
    sprintf('Exploding solutions found!')
    result = initconds([1,2,3,4],status == false);
    otherconds = initconds([1,2,3,4],status == true);
    figure(2)
    hold on
    scatter(result(1,:),result(2,:),'r','d')
    scatter(otherconds(1,:),otherconds(2,:),'b','*')
    xlabel('$y_1$','FontSize', 12, 'interpreter','latex') % x-axis label
    ylabel('$y_2$','FontSize', 12, 'interpreter','latex') % y-axis label
    legend('Unstable','Stable');
    hold off
    result
    condstested = [initconds(1,:);initconds(2,:)]';
end
