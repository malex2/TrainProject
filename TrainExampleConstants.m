close all
clear

M1 = 0;
M2 = 0;
k  = 4;
%F  = 1;
mu = 0.02;
g  = 9.8;

% Train PID 

t0 = 0; %Initial time of simulation in seconds
tf = 300; %Final time of simulation in seconds
h = 0.1; %time step

%Set gain
% kp = 0.06666666666667;
% ki = 0.01;
% kd = 0;

kp = 0.0;
ki = 0.3;
kd = 0.0;

% kp = 1;
% ki = 1;
% kd = 0;

% % Initialize Vectors
% n=(tf-t0)/h+1; %Size of time/all vectors
% error = zeros(1,n);
% Derror = zeros(1,n);
% Ierror = zeros(1,n);
% F = zeros(1,n);
% x1 = zeros(1,n);
% x2 = zeros(1,n);
% v1 = zeros(1,n);
% v2 = zeros(1,n);
% a1 = zeros(1,n);
% a2 = zeros(1,n);
% 
% % Set intial and desired values 
% %Desired velocity
% 
% % % Linear desired quantity
% % vd0 = 0; % For step input, vd0 = vd1
% % vd1 = 1;
% % t1 = 150;
% % y = @(x)( (vd1-vd0)/(t1-t0)*(x-t0)+vd0 );
% % vd = zeros(1,n);
% % vd(1:(t1-t0)/h+1) = y(1:h:t1+1);
% 
% % Step input of desired quantity
% Step = 1; % Magnitude of step
% t1 = 150; %Time of step down. t1 = tf for constant desired quantity
% vd = zeros(1,n);
% vd(1:(t1-t0)/h+1) = Step;
% 
% v1(1,1) = 0; %initial velocity of train
% v2(1,1) = 0; %initial velocity of traincar
% x1(1,1) = 0; %initial position of train
% x2(1,1) = 0; %initial position of traincar
% 
% i = 1;
% 
% for t = t0:h:tf
%   %Train car error
%   error(i) = vd(i) - v1(i);
%   
%   %Derivate of error
%    if i>1
%      Derror(i) = (error(i) - error(i-1))/h;
%    else
%      Derror(i) = 1;
%    end
%    
% %  if i>1
% %  Derror(i) = -(v1(i) - v1(i-1))/h;
% %  else
% %     Derror(i) = 1;
% %  end
% 
% if i>1
% Ierror(i) = Ierror(i-1) + error(i)*h;
% else
%     Ierror(i) = 0;
% end
%   
%   %Compute force
%    F(i) = kp*error(i) + kd*Derror(i) + ki*Ierror(i);
%    
%   %Evaluate System Dynaimcs
%   %1-train 2-train car
%   
%   a1(i+1) = 1/M1*(F(i) - (x1(i)-x2(i))*k - mu*g*M1*v1(i));
%   a2(i+1) = 1/M2*((x1(i)-x2(i))*k - mu*g*M2*v2(i));
%   
%   %Compute velocities (integral of a1 and a2)
%   
%   v1(i+1) = v1(i) + a1(i+1)*h;
%   v2(i+1) = v2(i) + a2(i+1)*h;
%   
%   %Compute positions (integral of v1 and v2)
% 
%   x1(i+1) = x1(i) + v1(i+1)*h;
%   x2(i+1) = x2(i) + v2(i+1)*h;
%   
%   i = i+1; 
% end
% 
% %Plot results
% 
% t = t0:h:tf;
% vd = vd.*ones(1,length(t));
% 
% figure
% plot(t,vd,' k',t,v1(1:end-1),' g')
% legend('Desired Velocity','Velocity')
% title('Train Velocity vs Time')
% ylabel('Velocity m/s')
% xlabel('Time s')
% axis([0 300 0 1.5])
% 
% % figure
% % plot(t,x1(1:end-1))
% % title('Train Position vs Time')
% % ylabel('Position m')
% % xlabel('Time s')
