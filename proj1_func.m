function [ coordinates ] = proj1_func(y1_0,y2_0,y3_0,y4_0,alpha,beta,nu,kappa,NumSteps )
% Uses Euler's method to solve our 4D differential equation, using the coordinate
% changed system. This function will output a graph and "coordinates" 
% vector for this trajectory's y_1 and y_2 coordinate only. As MATLAB
% cannot effectively display 4d space, we must choose which plane in 4d we
% wish to view. If we desire to view another plane, this function must be
% changed to view another pair of coordinates. 
% Finally, this function takes in a "kappa" value, for whether or not
% random noise is desired, with "1" indicating random noise, and "0"
% indicating no random noise. This function itself must be adjusted so that
% the random noise is in the correct direction (if kappa is set to 1).
% NumSteps is a measure of time we let it run for, where NumSteps = 1000
% would be a single unit of time.
% Takes in inputs (y1_0,y2_0,y3_0,y4_0,alpha,beta,nu,kappa,NumSteps )
% Example: proj1_func(2,0,0,2,1,1,1,0,2000 );
%

close all
rng('shuffle');
 
y1(1)=y1_0;%IC
y2(1)=y2_0;
y3(1)=y3_0;
y4(1)=y4_0;

% Estimate of dt necessary for stability of Euler's method
C = .001/max([log(y1(1)^2+y2(1)^2),1]); 
dt= .001*C; 
N= NumSteps/(1000*dt); %Calculates N to preserve total time (NumSteps*dt=5000)
sdt=sqrt(dt);


% Euler's method is run, and the resultant trajectory is recorded
for i=1:N
    db=sdt*randn;
    y1(i+1) = y1(i) + (-nu*y1(i)+beta*(y1(i)^2-y2(i)^2)-beta*(y3(i)^2-y4(i)^2))*dt;  
    y2(i+1) = y2(i) + (-nu*y2(i))*dt;
    y3(i+1) = y3(i) + (-nu*y3(i)+2*beta*(y1(i)*y3(i)-y2(i)*y4(i)))*dt + 1/2*db*kappa; %if kappa is set to 0, this term will have no effect
    y4(i+1) = y4(i) + (-nu*y4(i))*dt + 1/2*db*kappa;
end


% Displayed are the y_1, y_2 axes. If you would like to view the y_3, y_4
% coordinates on the same axes, simply uncomment the second plot statement
hold on
plot(y1,y2,'b',0,0,'.r',y1(1),y2(1),'.r', 'MarkerSize',20)
%plot(y3,y4,'r',0,0,'.r',y3(1),y4(1),'.r', 'MarkerSize',20)
xlabel('y1','FontSize', 12, 'interpreter','latex') % x-axis label
ylabel('y2','FontSize', 12, 'interpreter','latex') % y-axis label
str = sprintf('Trajectory of $y1$, $y2$ for alpha = %d, beta = %d, nu = %d',alpha,beta,nu);
title(str,'FontSize', 16, 'interpreter','latex')
coordinates = [y1;y2];
hold off
end
