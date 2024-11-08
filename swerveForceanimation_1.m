%Matlab code for swerve angle PID simulation 
%with field geometry and motor friction (amps free added as function of speed sign).
%coded by chris vamfun  mentor FRC team 599 10/24
% contact for questions vamfun@Yahoo.com
clear
clc
close all

%simulation time step dt and number of steps
dt =.001; % integration step size ..sec
T = 4  % run time ..sec

% command inputs to joy stick
vx_joy = 1, vy_joy = 0, omega_joy = 2*pi  % swerve velocity commands in m/s,omega commands in rad/s

fudge = 1.30 %used for experimental compensation 
%initialize pid error variables  

nstep =T/dt+1;
prnt_dt= T
pose_dt = .6% how often robot pose is plotted

swerve_cmd_loop_dt= .05;% set update rate for swerve cmd generation

loop_step= round(swerve_cmd_loop_dt/dt);
mstep = round(prnt_dt/dt);
posestep = round(pose_dt/dt); %How often show pose frame
iprnt=1;
iloop=1;
ipose=1; 
iframe = 1;

err_last=0; 
cumerr=0; 
%field state initialization
w=0*omega_joy;
psi =0 ;  % deg body heading angle , positive counterclockwise. 0 facing opposing alliance

Vx=0*vx_joy ; % x field speed fps
Vy=0 ; % y field speed fps 

X =0 ; %m positive going away from your alliance station 
Y= 0;  %m positive left

%robot constants here
rw = .0508 ;% 2 in radius of wheel 
% robot mass in kg
m = 70.31   % 155 lb

L=22.75/39.37; W = 22.75/39.37;xo=0 ;% xo is cg offset in x direction
wheel_positions = [
        xo+L/2, W/2;   % Front-left
        xo+L/2, -W/2;  % Front-right
        xo-L/2, -W/2;  % Back- right
        xo-L/2, W/2; % Back-left
    ];
%wheel cornering power.. cf_per_rad
cf_per_rad = m*9.8/4;% assumes coeff of friction = 1 and 25% weight on wheel

%robot spin moment of inertia kg m**2
%assume uniform density disc with radius of radius 14 in
% Izz = .5*rg^2*m     ,rg = .3556 m (14 in)
rg = .3556 
Izz= .5*rg*rg*m + xo^2*m;% Izz=4.445
%KRAKEN motor specs
motor_specs=[7.09, %stall torque nm
       6000,% free speed rpm
       366, % stall current amps
       2.0] % free current amps
   gear_v=6.75
   amp_limit = 40  % motor current limit
   %gear_alpha = 150/7
%swerve motor tau ,linear motion
tau= m*rw^2*motor_specs(2)*2*pi/60/motor_specs(1)/(gear_v)^2/4
tau_alpha= .05 %wheel angle response time constant
%exp first order filter constants c1,c2
%used to model wheel angle servo lags
c1 = exp(-dt/tau_alpha);
c2 = 1 - c1;
% z = z*c1 + cmd*c2 ,

%swerve wheel velocity control PID constants 
kff=1
Kp =20 %nominal value 10 
Kd = 0%tau*Kp % tau = .054 s
amps=0.;%motor current amps
nm=0;%motor torq nm
%swerve outer loop velocity control PI constants
Kp_swerve =0;
Ki_swerve = 0;
delay_dt_sec = fudge* swerve_cmd_loop_dt;
cumerr_vx=0; cumerr_vy = 0;
% maximum wheel speed
vwmax=motor_specs(2)*2*pi/60*rw/gear_v
% predivide Kd by dt 
Kd=Kd/dt; 
%output array
time=zeros(1,nstep);
out1 =zeros(1,nstep);
out2 =zeros(1,nstep);
out3 =zeros(1,nstep);
out4 =zeros(1,nstep);
out5 =zeros(1,nstep);
out6 =zeros(1,nstep);
out7 =zeros(1,nstep);
out8 =zeros(1,nstep);
out9 =zeros(1,nstep);
out10 =zeros(1,nstep);
outX = zeros(1,nstep);
outY = zeros(1,nstep);
outPsi= zeros(1,nstep);
alphaw = zeros(1,4);
alpha = zeros(1,4);
beta = zeros(1,4);
alphac = zeros(1,4);
t = 0; % numerical integration loop *****************
figure
for i = 1 :nstep
   if i == iloop  % update command loop 
        iloop = iloop + loop_step;
%generate wheel speed and angle commands from joy stick. also speed and angle rates for kicks

% outer loop field velocity PID control 
       
         delta =-w*delay_dt_sec/2; % drift angle (rad) correction for time delays
         joy_vx_ff= vx_joy - vy_joy*sin(delta);
         joy_vy_ff= vy_joy + vx_joy*sin(delta);
         err_vx= (vx_joy- Vx);
         err_vy= (vy_joy- Vy);
         cumerr_vx = cumerr_vx + err_vx*swerve_cmd_loop_dt;
         cumerr_vy = cumerr_vy + err_vy*swerve_cmd_loop_dt;
      
         vx_pi = joy_vx_ff + Kp_swerve*(err_vx) + Ki_swerve*cumerr_vx;
         vy_pi = joy_vy_ff + Kp_swerve*(err_vy) + Ki_swerve*cumerr_vy;
       

       [vwc, alphac ,v_wheels_ratec,alpha_wheels_ratec] = ComputeSwerveCommandsfromJoy(vx_pi, vy_pi,omega_joy, psi, alphac,wheel_positions,vwmax);

   end
      
          
      [vw, alpha] = fieldOrientedToWheelVelocities(Vx, Vy, w, psi , alpha,wheel_positions);
          
             %wheel speed loop PID  calculations ********************
    err=1*(vwc-vw)/vwmax ; 
    errd = err-err_last;   
    err_last=err;
   
    cumerr=cumerr + err*dt;
    cmd= err*Kp + errd*Kd + kff*vwc/vwmax; %
    
 %wheel angle PID servo output simulation..simple lag
 %isnt necessary to code pid for angle loop
        alphaw= alphaw*c1 + alphac*c2;
    
    %sum forces and moments

   
        % wheel velocity motor torques and currents 
         
        for j = 1:4
        [nm(j),amps(j)]=MOTOR(cmd(j),vw(j)/rw,amp_limit,gear_v,1,motor_specs);
         end
        
        fw = nm/rw;%wheel forces
        beta = atan2(sin(alpha-alphaw),cos(alpha-alphaw))  ;% compute wheel side slip angle
        fwn= -cf_per_rad*sin(beta);% side force due to sideslip
        %rx and ry are column vectors :  fw,fwn,fx,fy are row vectors
        rx = wheel_positions(:,1)*cos(psi) - wheel_positions(:,2)*sin(psi);
        ry = wheel_positions(:,1)*sin(psi) + wheel_positions(:, 2)*cos(psi);
        fx = fw.*cos(alphaw + psi) -fwn.*sin(alphaw + psi);
        fy = fw.*sin(alphaw+ psi) +fwn.*cos(alphaw+ psi);
        
        %forces and moments

        
        moment= sum(-ry'.*fx + rx'.*fy);
        fxnet = sum(fx);
        fynet = sum(fy);
    
    wd = moment/Izz;
    ax = fxnet/m;
    ay = fynet/m;
    % odometry 
    if i > 1
    w = w + wd*dt;
    Vx= Vx + ax*dt;
    Vy= Vy + ay*dt;
    psi = psi + w*dt;
    X= X + Vx*dt;
    Y= Y + Vy*dt;
    t = t + dt ;
    end
    time(i)= t;
    out10(i) = alphaw(4);
    out9(i) = alphaw(3);
    out8(i) = alphaw(2);
    out7(i) = alphaw(1);
    out6(i)= beta(1);
    out5(i) = alphaw(1);
    out4(i) = vw(1);
    out3(i) = err(1); 
    out2(i) = err_vy;
    out1(i) =vw(1)/vwmax; 
    outX(i) = X;
    outY(i) =Y;
    outPsi(i) = psi;
    if i == iprnt
        iprnt = iprnt + mstep;
        
        %t,vwc,vw,w,Vx,Vy ,psi
        %alphaw+psi ,beta 
        %v_wheels_ratec,alphac
        %t,alphac ,alphaw,alpha,beta
      t,X,Y
    end
  if i == ipose
      %clf;
        ipose = ipose + posestep;  
       
    plotRobotpose(outX(i),outY(i),wheel_positions,outPsi(i),vw/vwmax,alphaw,.5)
    iframe = iframe + 1
    movieVector(iframe) = getframe;
    pause(.5)
    xlabel('x (m)');
ylabel('y (m)');
title('Trajectory of the Swerve Drive Robot ','with command compensation');
grid on;
axis equal;

  end

end  




plot(outX, outY, 'b-', 'LineWidth', 2);


%hold on

% xlabel('x (m)');
% ylabel('y (m)');
% title('Trajectory of the Swerve Drive Robot');
% grid on;
% axis equal;

if 0 % set to one to plot

% Plot the velocities (optional, if supported)
figure;
time = linspace(0, T, nstep);
%subplot(2, 1, 1);
plot(time, out2, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Vy field velocity error m/s ');
title('Vy field velocity err over Time');
grid on;
end
if 0 % set to one to plot

% Plot the velocities (optional, if supported)
figure;
time = linspace(0, T, nstep);
subplot(2, 1, 1);
plot(time, out3, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('wheel velocity error 1 m/s ');
title('wheel velocity err over Time');
grid on;

subplot(2, 1, 2);
plot(time, out4, 'g-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('wheel velocity (m/s)');
title('wheel velocity over Time');
grid on;
end
if 0 % set to 1 to plot

% Plot alpha and beta(optional, if supported)
figure;
time = linspace(0, T, nstep);
subplot(2, 1, 1);
plot(time, out5, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('alpha 1 (rad)');
title('wheel angle  over Time');
grid on;

subplot(2, 1, 2);
plot(time, out6, 'g-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('beta 1 (rad)');
title('wheel slip angle over Time');
grid on;
end

if 0 %set to 1 to plot

% Plot alpha wheel (optional, if supported)
figure;
time = linspace(0, T, nstep);
subplot(4, 1, 1);
plot(time, out7, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('alpha 1 (rad)');
title('wheel angle  over Time');
grid on;

subplot(4, 1, 2);
plot(time, out8, 'g-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('alpha 2 (rad)');
title('wheel angle over Time');
grid on;

time = linspace(0, T, nstep);
subplot(4, 1, 3);
plot(time, out9, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('alpha 3 (rad)');
title('wheel angle  over Time');
grid on;

subplot(4, 1, 4);
plot(time, out10, 'g-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('alpha 4 (rad)');
title('wheel angle over Time');
grid on;
end


function [deg360 deg180]= angle360or180(deg)
        if deg > 180 
            deg180 = deg-360;
        elseif deg < -180 
            deg180 = deg + 360;
        else 
            deg180 = deg;
        end
				if deg180 < 0
    				deg360 = deg180 + 360;
				else 
    				deg360 = deg180;
				end
end 
       
function [deg360, deg180]= ATAN2D(x,y)
  
        if x<0 
            deg360 = 180 - atand(-y/x);
        elseif(x>0) 
            if(y>0) 
                deg360 = atand(y/x);
            else 
                deg360 = 360 + atand(y/x);
            end
        elseif (y>0 ) 
            deg360= 90;
        else 
            deg360=270;
        end
       if deg360 >180 
           deg180 = deg360 -360;
       else deg180= deg360;
       end
end
function [nm,ampsl] = MOTOR(cmd,w_rps,amp_limit,gr,num,specs)  % w_rps is for wheel , i.e. gear box output rotation rate.  torque output is after gearbox
eff = .95;  %gear train efficiency
nm_stall= specs(1);
 amp_stall=specs(3) ;
rpm_free= specs(2);
rps_free= rpm_free*2*pi/60;
amp_free= specs(4);
amp_frict= 3 ;%  current required to overcome system rotational friction (e.g. gearbox)
ampsl=0;
ampslf=0;

if cmd>1 
    cmdl = 1;
elseif cmd< -1 
    cmdl= -1 ;
else cmdl=cmd;
end

amps=amp_stall*(cmdl- w_rps*gr/rps_free) ;%these amps don't include amps_free when @ Rpm_free
if amps>amp_limit     
    ampsl=amp_limit;
elseif amps<-amp_limit 
    ampsl=-amp_limit;
else 
    ampsl= amps ;
end
% add friction torque losses by subtracting equivalent current bias  amp_frict ; always opposes direction of rotation.
if w_rps>0
    ampslf=ampsl-amp_frict;
elseif w_rps < 0 
    ampslf= ampsl+amp_frict;
else ampslf = ampsl;
end
    nm= eff*num*gr*nm_stall*ampslf/(amp_stall-amp_free);


end

    
function [v_wheels, alpha_wheels ] = fieldOrientedToWheelVelocities(vx_field, vy_field, omega, theta,alpha_wheels_last ,wheel_positions)
    %convert field velocities to body axes
    vx_body = vx_field*cos(theta) + vy_field*sin(theta);
    vy_body = -vx_field*sin(theta) + vy_field*cos(theta);
   
    % front of robot is pointed along positive x axis   CCW rotation is +
     
  
   for i = 1:4
    % Calculate wheel position relative to the robot center   
    pos_x = wheel_positions(i, 1);
    pos_y = wheel_positions(i, 2);
  
   
    vx = vx_body - omega * pos_y ;
    vy = vy_body + omega * pos_x ;
    v_wheels(i) = sqrt(vx^2 + vy^2); 
    % Calculate wheel steering angles
    alpha_wheels(i) = atan2(vy, vx); %
 
   
   end
   %alpha_wheels = continuousAngle(alpha_wheels,alpha_wheels_last);
 end 
    function [v_wheelsc, alpha_wheelsc ,v_wheels_ratec,alpha_wheels_ratec] = ComputeSwerveCommandsfromJoy(vx_joy, vy_joy,omega_joy ,theta, alpha_wheels_last,wheel_positions,vwmax)
         max_temp=vwmax;
           %convert field velocities to body axes
    vx_body = vx_joy*cos(theta) + vy_joy*sin(theta);
    vy_body = -vx_joy*sin(theta) + vy_joy*cos(theta);
   
    % front of robot is pointed along positive x axis   CCW rotation is +
     
  
   for i = 1:4
    % i th wheel position relative to the robot center   
    pos_x = wheel_positions(i, 1);
    pos_y = wheel_positions(i, 2);
  
   
    vx = vx_body - omega_joy * pos_y ;
    vy = vy_body + omega_joy * pos_x ;
    v_wheels(i) = sqrt(vx^2 + vy^2); 
    % Calculate wheel steering angles
    alpha_wheelsc(i) = atan2(vy, vx); %
          % Calculate wheel accelerations and angle rates in
    % body axes caused by omega assumining vx_field, vy_field and omega constant
    ax =- pos_x*omega_joy^2;
    ay =- pos_y*omega_joy^2;  
    v_wheels_ratec(i) = (ax*vx + ay*vy)/v_wheels(i);
    alpha_wheels_ratec(i) = (-ax*vy + ay*vx)/v_wheels(i)^2;

  if v_wheels(i) > max_temp
        max_temp = v_wheels(i);
  end
  end
   % if overspeed then scale all speeds by factor that reduces highest speed to = vwmax.
   if max_temp > vwmax 
        factor = vwmax/max_temp;
        v_wheelsc=v_wheels*factor;
        %v_wheels_rate=v_wheels_rate/factor;
   else
       v_wheelsc = v_wheels;
    end 
    
       %alpha_wheelsc = continuousAngle(alpha_wheelsc,alpha_wheels_last);
   
end

function rate = smoothrate(a,alast,dt)
    % corrects for angle going from pi to -pi or -pi to pi
    % before deriving rate.
    delta = [0,0,0,0];
    delta = a - alast;
    n= size(a,2);
    for i = 1 : n
    if delta(1,i) > pi
        delta(1,i) = delta(1,i) - 2*pi;
   elseif delta(1,i) < -pi
        delta (1,i) = delta(1,i)  + 2*pi;
    end
    end 
    rate = delta/dt;
    end
function angle = continuousAngle(a,alast)
    % corrects for angle going from pi to -pi or -pi to pi .... makes angle continous
    
    delta = a - alast;
    m=round(delta/pi);
  
   angle = a - m*pi;
    
end

function plotRobotpose(X,Y,wheel_positions,theta,norm_vel,alphas,mag)
    [xp,yp]=createVerticies(X,Y,wheel_positions,theta);
    plot(xp,yp,"m");
grid on;
axis equal;

hold ("on")
   [velx,vely] = plotVelocityVectors(xp,yp,norm_vel,alphas,theta,mag);
   % plot robot index line
   xindex = [X + wheel_positions(1,1)*cos(theta)];
   yindex = [Y + wheel_positions(1,1)*sin(theta)];
   plot(xindex,yindex,"mo","LineWidth",2)
end

function [xp,yp]= createVerticies(x0,y0,points,theta)
    xp = zeros(1,5);
    yp = zeros(1,5);
    %rotate body to field 
    R = [ cos(theta) , -sin(theta) 
          sin(theta) , cos(theta)];
        
          xyrot = R*[points(:,1),points(:,2)]';
    p= [x0,y0] + xyrot';
    xp=[p(1,1), p(2,1) ,p(3,1) ,p(4,1) ,p(1,1)];
    yp=[p(1,2), p(2,2) ,p(3,2) ,p(4,2) ,p(1,2)];
    
end 

 function [velx,vely] = plotVelocityVectors(xp,yp,norm_vel,alphas,theta,mag)
     %this function attaches a velocity vector to each wheel vertex
     % mag parameter scales the length of the velocity vector arrow 
for i = 1:4
    xv(i) = xp(i) + norm_vel(i)*mag*cos(alphas(i)+theta);
    yv(i)  = yp(i) + norm_vel(i)*mag*sin(alphas(i)+theta);
    velx = [xp(i) , xv(i)];
    vely = [yp(i) , yv(i)];
    plot(velx,vely,"g","LineWidth", 2);
end
 end
    
         
         
 