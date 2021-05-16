clear, clc, clf

%{
    z (up)
    |
    |
    |_______ y (towards right)
   / 
  / 
 x (towards back)

vehicle ground plane at z = 0

Roll center not working?
- read roll center myths and reality, computing suspension books, RCVD

TODO in next 2 weeks:
- verify current outputs with solidworks
- check roll center and roll calculations
- document how roll is found (approximation)
- anti-stuff (dive/squat)
LATER:
- bonus lateral roll center position
- % ackermann for steering geometry
- steering wheel to tire sensitivity 
%}

% wheelcenter
wc_orig     = [522.8;   614.05;	190.6]; 

% outboard points (ob)
upper       = [522.8;	517.2;	262.75];
lower       = [450.2;	529.9;	123.59];
tie_ob      = [570.0;	529.9;	115.61];

% inboard points (ib)
upper_f     = [336.0;	280.7;	249.0];
upper_r     = [565.0;	280.8;	240.6];
lower_f     = [274.4;	260.4;	136.3];
lower_r     = [567.0;	167.0;	117.0];
tie_ib      = [591.2;	150.3;	115.5];

% initial values
init_camber_deg  = -1.5;
init_toe_deg     = 0;
init_track_mm    = wc_orig(2);

% iteration constants
loops            = 100;    % must be even - half upwards, half downwards
time_interval    = 0.5;    % step size in seconds 
input_vel        = 1;      % upwards wc velocity in mm/s
filename = 'vector_movement.gif';

% compiling starting points and saving original values
wc          = wc_orig;
wc3_orig    = [wc, wc, wc];
wc3         = wc3_orig;
ob3_orig    = [upper, lower, tie_ob];
ob3         = ob3_orig;
ob5_orig    = [ob3(:,1), ob3(:,1), ob3(:,2), ob3(:,2), ob3(:,3)];
ob5         = ob5_orig;
ib5_orig    = [upper_f, upper_r, lower_f, lower_r, tie_ib];
Rx_          = [1 0 0; 0 cos(init_camber_deg*pi/180) -sin(init_camber_deg*pi/180);...
             0 sin(init_camber_deg*pi/180) cos(init_camber_deg*pi/180)];
contact_orig      = Rx_*[0; 0; -1 * wc_orig(3)];
wc_to_contact     = contact_orig;

% set up data collection arrays
vertical_data   = zeros(1, loops+1);
wheelbase_data  = zeros(1, loops+1);
htrack_data     = zeros(1, loops+1); % half of track (single wheel)
camber_data     = zeros(1, loops+1);
spin_data       = zeros(1, loops+1);
toe_data        = zeros(1, loops+1);
htrack_data(1)  = init_track_mm;
camber_data(1)  = init_camber_deg;
toe_data(1)     = init_toe_deg;
roll_data       = zeros(1, loops+1);
rollheight_data = zeros(1, loops+1);

cp_to_ic        = zeros(3, loops+1);
cont_data       = zeros(3, loops+1);
ca_data         = zeros(15,loops+1);
ob_to_wc_data   = zeros(3, loops+1);

wc_to_ic = zeros(3, loops+1);
ctp_to_ic = zeros(3, loops+1);

% set up variables to track change
vertical_mm     = 0;
wheelbase_mm    = 0;
track_mm        = init_track_mm; 
camber_deg      = init_camber_deg;
spin_deg        = 0;
toe_deg         = init_toe_deg;
roll_deg        = 0;

% find plotting limits
max_x = max([ib5_orig(1,:), ob5_orig(1,:), wc(1)]) + 10;
min_x = min([ib5_orig(1,:), ob5_orig(1,:), wc(1)]) - 10;
max_y = max([ib5_orig(2,:), ob5_orig(2,:), wc(2)]) + 10;
min_y = min([ib5_orig(2,:), ob5_orig(2,:), wc(2)]) - 10;
max_z = max([ib5_orig(3,:), ob5_orig(3,:), wc(3)]) + 50;
min_z = min([ib5_orig(3,:), ob5_orig(3,:), wc(3)]) - 50;

% suppress output during runtime
%fig = figure;
%fig.Visible = 'off';

% set up gif image array
mov(loops) = struct('cdata', [], 'colormap', []);
mov2(loops+1) = struct('cdata', [], 'colormap', []);

 for iter = 1:loops    
    %    Vx         Vy     Wx      Wy    Wz
    syms dWheelbase dTrack dCamber dSpin dToe
    
    % linkage vectors: 1=FUCA, 2=RUCA, 3=FLCA, 4=RLCA, 5=TROD
    ca5 = ib5_orig - ob5;
    
    % outputting control arm lengths
    fprintf('ca1: %f\nca3: %f\nca5: %f\n\n', ...
        norm(ca5(:,1)), norm(ca5(:,3)), norm(ca5(:,5)));
    
    % radius from wheelcenter
    rad3 = ob3 - wc3;
    
    % velocity vector at wheelcenter (x3)
    wc_veloc = [dWheelbase, dWheelbase, dWheelbase; 
                dTrack,     dTrack,     dTrack;...
                input_vel,  input_vel,  input_vel];
    
    % relative velocities at outboard points
    wc_angular = [dCamber; dSpin; dToe];
    wc_angular3 = [wc_angular, wc_angular, wc_angular];
    rel_veloc = cross(wc_angular3, rad3);
    
    % velocities at outboard points: 1=Upper, 2=Lower, 3=Tie
    ob_veloc = rel_veloc + wc_veloc;
    
    % dot product of velocities with each control arm
    eqn1 = dot(ob_veloc(:,1), ca5(:,1)) == 0;
    eqn2 = dot(ob_veloc(:,1), ca5(:,2)) == 0;
    eqn3 = dot(ob_veloc(:,2), ca5(:,3)) == 0;
    eqn4 = dot(ob_veloc(:,2), ca5(:,4)) == 0;
    eqn5 = dot(ob_veloc(:,3), ca5(:,5)) == 0;
    
    % solve linear equation
    sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5],...
        [dWheelbase, dTrack, dCamber, dSpin, dToe]);
    
    % subsitute variables
    ob_veloc_sub = subs(ob_veloc,...
                        [dWheelbase,     dTrack,     dCamber,     dSpin,     dToe], ...
                        [sol.dWheelbase, sol.dTrack, sol.dCamber, sol.dSpin, sol.dToe]);
    wc_veloc_sub = subs(wc_veloc(:,1),...
                        [dWheelbase,     dTrack],...
                        [sol.dWheelbase, sol.dTrack]);
    
    % integrate wheelcenter velocities
    dt = time_interval;
    % velocity in coordinate unit(mm) per second --> mm
    ob_displ = ob_veloc_sub*dt;
    wc_displ = wc_veloc_sub*dt;
    % angular velocity in rad per second --> rad
    ang_vel = [sol.dCamber; sol.dSpin; sol.dToe];
    ang_displ = ang_vel*dt;
    
    % record data
    wheelbase_mm    = wheelbase_mm + wc_displ(1);
    track_mm        = track_mm + wc_displ(2);
    vertical_mm     = vertical_mm + wc_displ(3);
    camber_deg      = camber_deg + -1 * ang_displ(1)*180/pi;
    spin_deg        = spin_deg + ang_displ(2)*180/pi;
    toe_deg         = toe_deg + ang_displ(3)*180/pi;
    
    vertical_data(iter+1)   = vertical_mm;
    wheelbase_data(iter+1)  = wheelbase_mm;
    htrack_data(iter+1)     = track_mm;
    camber_data(iter+1)     = camber_deg;
    spin_data(iter+1)       = spin_deg;
    toe_data(iter+1)        = toe_deg;
    
    % update points
    Rx = [1 0 0; 0 cos(ang_displ(1)) -sin(ang_displ(1)); 0 sin(ang_displ(1)) cos(ang_displ(1))];
    Ry = [cos(ang_displ(2)) 0 sin(ang_displ(2)); 0 1 0; -sin(ang_displ(2)) 0 cos(ang_displ(2))];
    Rz = [cos(ang_displ(3)) -sin(ang_displ(3)) 0; sin(ang_displ(3)) cos(ang_displ(3)) 0; 0 0 1];
    
    rad3 = Rx*rad3;
    rad3 = Ry*rad3;
    rad3 = Rz*rad3;
    
    wc_to_contact = vpa(Rx*wc_to_contact);
    vel_contact = cross(ang_vel, wc_to_contact) + wc_veloc_sub;
    syms v_dist
    sol2 = solve(dot([0; -1; v_dist], vel_contact) == 0, v_dist);
    cp_to_ic(:,iter+1) = [0; -1; sol2]./norm([0; -1; sol2]);
    
    wc_to_ic(:,iter+1) = [1 0 0; 0 cos(pi/2) -sin(pi/2); 0 sin(pi/2) cos(pi/2)]*[0;wc_veloc_sub(2);wc_veloc_sub(3)];
    ctp_to_ic(:,iter+1) = [1 0 0; 0 cos(pi/2) -sin(pi/2); 0 sin(pi/2) cos(pi/2)]*[0;vel_contact(2);vel_contact(3)];
    cont_data(:,iter+1)     = wc_to_contact;
    
    % plot vectors
    quiver3(ob5(1,:), ob5(2,:), ob5(3,:), ca5(1,:), ca5(2,:), ca5(3,:), 'AutoScale', 'off');
    %view(90, 0)
    hold on
    quiver3(wc3(1,:), wc3(2,:), wc3(3,:), rad3(1,:), rad3(2,:), rad3(3,:), 'AutoScale', 'off');
    quiver3(wc(1), wc(2), wc(3), wc_to_contact(1), wc_to_contact(2), wc_to_contact(3), 'AutoScale', 'off');
    contact = wc + wc_to_contact;
    quiver3(contact(1), contact(2), contact(3), vel_contact(1)*50, vel_contact(2)*50, vel_contact(3)*50, 'AutoScale', 'off');
    quiver3(contact(1), contact(2), contact(3), cp_to_ic(1,iter+1)*100, cp_to_ic(2,iter+1)*100, cp_to_ic(3,iter+1)*100, 'AutoScale', 'off');
    hold off
    set(gca, 'DataAspectRatio', [1 1 1],...
        'XLim', [min_x max_x],...
        'YLim', [-50 max_y],...
        'ZLim', [-50 max_z])
    drawnow limitrate
    mov(iter) = getframe;
    
    % make a GIF file
    im = frame2im(mov(iter));
    [imind, cmap] = rgb2ind(im,256);
    if iter == 1
        imwrite(imind, cmap, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
    elseif mod(iter, 5) == 0
        imwrite(imind, cmap, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
    
    % only move wheel center and rotate radius vectors
    wc  = vpa(wc + wc_displ); % vpa( ) --> round for efficiency
    wc3 = [wc, wc, wc];
    ob3 = vpa(wc3 + rad3);
    ob5 = [ob3(:,1), ob3(:,1), ob3(:,2), ob3(:,2), ob3(:,3)];
    
    ca5 = ib5_orig - ob5;
    ca_data(1:3,iter+1)     = -ca5(:,1);
    ca_data(4:6,iter+1)     = -ca5(:,2);
    ca_data(7:9,iter+1)     = -ca5(:,3);
    ca_data(10:12,iter+1)   = -ca5(:,4);
    ca_data(13:15,iter+1)   = -ca5(:,5);
    ob_to_wc_data(:,iter+1) = -rad3(:,1);
    
    % flip from upwards input velocity to downwards
    if iter == loops/2
        input_vel = -1 * input_vel;
        wc     = wc_orig;
        wc3    = wc3_orig;
        ob5    = ob5_orig;
        ob3    = ob3_orig;
        vertical_mm    = 0;
        wheelbase_mm   = 0;
        track_mm       = init_track_mm;
        camber_deg     = init_camber_deg;
        spin_deg       = 0;
        toe_deg        = init_toe_deg;
        roll_deg       = 0;
        wc_to_contact        = contact_orig;
    end
 end

%fig.Visible = 'on';
%movie(mov);

vertical_data_ordered   = [flip(vertical_data(loops/2+2:loops+1)),...
                            vertical_data(1:loops/2+1)];
htrack_data_ordered     = [flip(htrack_data(loops/2+2:loops+1)),...
                            htrack_data(1:loops/2+1)];
camber_data_ordered     = [flip(camber_data(loops/2+2:loops+1)),...
                            camber_data(1:loops/2+1)];
toe_data_ordered        = [flip(toe_data(loops/2+2:loops+1)),...
                            toe_data(1:loops/2+1)];
roll_data_ordered       = [flip(roll_data(loops/2+2:loops+1)),...
                            roll_data(1:loops/2+1)];             
cp_to_ic_right          = [flip(cp_to_ic(:,loops/2+2:loops+1), 2),...
                            cp_to_ic(:,1:loops/2+1)];
% setting initial roll center to be first iteration
cp_to_ic_right(:,51)    = cp_to_ic_right(:,50);
% invert order for left side  
cp_to_ic_left           = flip(cp_to_ic_right, 2);
% flip y direction to left side
cp_to_ic_left(2,:)      = -1*cp_to_ic_left(2,:);

wc_to_ic                =[flip(wc_to_ic(:,loops/2+2:loops+1), 2),...
                            wc_to_ic(:,1:loops/2+1)];
ctp_to_ic                =[flip(ctp_to_ic(:,loops/2+2:loops+1), 2),...
                            ctp_to_ic(:,1:loops/2+1)];
                        
ca_ordered         = [flip(ca_data(:,loops/2+2:loops+1), 2),...
                            ca_data(:,1:loops/2+1)];
ob_to_wc_ordered   = [flip(ob_to_wc_data(:,loops/2+2:loops+1), 2),...
                            ob_to_wc_data(:,1:loops/2+1)];
cont_ordered       = [flip(cont_data(:,loops/2+2:loops+1), 2),...
                            cont_data(:,1:loops/2+1)];
figure2 = figure;
filename2 = 'roll_center.gif';

cont_ordered(:,loops/2+1) = contact_orig;
ob_to_wc_ordered(:,loops/2+1) = ob_to_wc_ordered(:,loops/2+2);
ca5_orig = ob5_orig - ib5_orig;
ca_ordered(:,loops/2+1) = [ca5_orig(:,1); ca5_orig(:,2); ca5_orig(:,3); ca5_orig(:,4); ca5_orig(:,5)];
init_ground = 0;
ib5 = ib5_orig;

% ROLL CALCULATIONS IN DEVELOPMENT:
for iter = loops/2+1:loops+1
    % find initial roll center
    syms rscale lscale
    rightV = [0; htrack_data_ordered(iter); 0] + (cp_to_ic_right(:,iter) .* rscale);
    leftV  = [0; -1*htrack_data_ordered(iter); 0] + (cp_to_ic_left(:,iter) .* lscale);
    rceqn2 = rightV(2) - leftV(2) == 0;
    rceqn3 = rightV(3) - leftV(3) == 0;
    rc_sol = solve([rceqn2, rceqn3], [rscale, lscale]);
    rc_vec = vpa(subs(rightV, rscale, rc_sol.rscale));
    rollheight_data(iter) = rc_vec(3);
    rightVisual = vpa(cp_to_ic_right(:,iter) .* rc_sol.rscale);
    leftVisual  = vpa(cp_to_ic_left(:,iter) .* rc_sol.lscale);
    rc_vec5 = [rc_vec, rc_vec, rc_vec, rc_vec, rc_vec];
    
    % rotate ob5 around roll center (increment system then make contact patch touch ground)
    syms theta
    Rx_ = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    rc_to_ib = ib5 - rc_vec5;
    rc_to_cp = rc_to_ib(:,1) + ca_ordered(1:3,iter) + ob_to_wc_ordered(:,iter) + cont_ordered(:,iter);
    rc_to_cp = Rx_*rc_to_cp;
    cp = rc_vec + rc_to_cp;
    new_sol = 0;
    if iter == loops/2+1
        cp = subs(cp, theta, 0);
        rc_to_cp = subs(rc_to_cp, theta, 0);
        Rx_ = subs(Rx_, theta, 0);
        init_ground = cp(3);
        disp(init_ground);
        roll_data_ordered(iter) = 0;
    else
        rotate_eqn = cp(3) == init_ground;
        rotate_sol = solve(rotate_eqn, theta);
        [new_sol, ind] = min([abs(rotate_sol(1)), abs(rotate_sol(2))]);
        new_sol = real(rotate_sol(ind)/abs(rotate_sol(ind))*new_sol);
        rc_to_cp = subs(rc_to_cp, theta, new_sol);
        Rx_ = subs(Rx_, theta, new_sol);
        cp = rc_vec + rc_to_cp;
        roll_data_ordered(iter) = roll_data_ordered(iter-1) + vpa(new_sol)*180/pi;
    end
    
    %update points
    ib5 = rc_vec5 + Rx_*rc_to_ib;
    
    % rotate iterations
    Rx_not = [1 0 0; 0 cos(-new_sol) -sin(-new_sol); 0 sin(-new_sol) cos(-new_sol)];
    for count = iter:loops+1
        ca_ordered(1:3, count)   = Rx_*ca_ordered(1:3, count);
        ca_ordered(4:6, count)   = Rx_*ca_ordered(4:6, count);
        ca_ordered(7:9, count)   = Rx_*ca_ordered(7:9, count);
        ca_ordered(10:12, count) = Rx_*ca_ordered(10:12, count);
        ca_ordered(13:15, count) = Rx_*ca_ordered(13:15, count);
        cont_ordered(:,count)    = Rx_*cont_ordered(:,count);
        ob_to_wc_ordered(:,count)= Rx_*ob_to_wc_ordered(:,count);
        wc_to_ic(:,count) = Rx_*wc_to_ic(:,count);
        ctp_to_ic(:,count) = Rx_*ctp_to_ic(:,count);
        if count ~= loops+1
            cp_to_ic_right(:,count+1)  = Rx_*cp_to_ic_right(:,count+1);
            cp_to_ic_left(:,count+1)   = Rx_*cp_to_ic_left(:,count+1); % <--- NOT ROLLED AROUND RC??
        end
    end
    
    
    ca_new5 = [ca_ordered(1:3,iter), ca_ordered(4:6,iter), ca_ordered(7:9,iter), ca_ordered(10:12,iter), ca_ordered(13:15,iter)]; 
    % visual output
    quiver3(rc_vec5(1,:), rc_vec5(2,:), rc_vec5(3,:), rc_to_ib(1,:), rc_to_ib(2,:), rc_to_ib(3,:), 'AutoScale', 'off', 'Color', 'y');
    view(90, 0)
    hold on
    quiver3(0,0,0, rc_vec(1), rc_vec(2), rc_vec(3), 'AutoScale', 'off', 'Color', 'g');
    quiver3(rc_vec5(1,:) + rc_to_ib(1,:), rc_vec5(2,:) + rc_to_ib(2,:), rc_vec5(3,:) + rc_to_ib(3,:), ...
        ca_new5(1,:), ca_new5(2,:), ca_new5(3,:), 'AutoScale', 'off', 'Color', 'b');
    quiver3(rc_vec(1) + rc_to_ib(1,1) + ca_new5(1,1) + ob_to_wc_ordered(1,iter),...
            rc_vec(2) + rc_to_ib(2,1) + ca_new5(2,1) + ob_to_wc_ordered(2,iter),...
            rc_vec(3) + rc_to_ib(3,1) + ca_new5(3,1) + ob_to_wc_ordered(3,iter), ...
            1000*wc_to_ic(1,iter), 1000*wc_to_ic(2,iter), 1000*wc_to_ic(3,iter), 'AutoScale', 'off', 'Color', 'm');
   quiver3(cp(1), cp(2), cp(3), 1000*ctp_to_ic(1,iter), 1000*ctp_to_ic(2,iter), 1000*ctp_to_ic(3,iter), 'AutoScale', 'off', 'Color', 'm');
        
    quiver3(rc_vec(1), rc_vec(2), rc_vec(3), rc_to_cp(1), rc_to_cp(2), rc_to_cp(3), 'AutoScale', 'off', 'Color', 'r');
    %quiver3(cp2(1), cp2(2), cp2(3), rightVisual(1), rightVisual(2), rightVisual(3), 'AutoScale', 'off');
    %quiver3(0, htrack_data_ordered(iter), 0, 100*rightVisual(1), 100*rightVisual(2), 100*rightVisual(3), 'AutoScale', 'off');
    %quiver3(0, -1*htrack_data_ordered(iter), 0, leftVisual(1), leftVisual(2), leftVisual(3), 'AutoScale', 'off');
    hold off
    set(gca, 'DataAspectRatio', [1 1 1],...
        'XLim', [0 1000],...
        'YLim', [-1250 700],...
        'ZLim', [-100 500])
    drawnow limitrate
    mov2(iter) = getframe;
    % make a GIF file
    im = frame2im(mov2(iter));
    [imind, cmap] = rgb2ind(im,256);
    if iter == loops/2+1
        imwrite(imind, cmap, filename2, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
    elseif mod(iter, 3) == 0
        imwrite(imind, cmap, filename2, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

roll_data_ordered(1:loops/2) = -flip(roll_data_ordered(loops/2+2:loops+1), 2);
rollheight_data(1:loops/2) = flip(rollheight_data(loops/2+2:loops+1), 2);
rollheight_data(loops/2+1) = rollheight_data(loops/2);

figure('Name', 'Camber vs Wheel Travel');
plot(vertical_data_ordered, camber_data_ordered, 'r');
title('Camber vs Wheel Travel')
xlabel('Wheel Travel (mm)')
ylabel('Camber Angle (degrees)')

figure('Name', 'Toe vs Wheel Travel')
plot(vertical_data_ordered, toe_data_ordered, 'r');
title('Toe vs Wheel Travel')
xlabel('Wheel Travel (mm)')
ylabel('Toe Angle (degrees)')

figure('Name', 'Camber vs Roll');
plot(roll_data_ordered, camber_data_ordered + roll_data_ordered, 'r');
title('Camber vs Roll')
xlabel('Roll Angle (degrees)')
ylabel('Camber Angle (degrees)')

figure('Name', 'Toe vs Roll');
plot(roll_data_ordered, toe_data_ordered, 'r');
title('Toe vs Roll')
xlabel('Roll Angle (degrees)')
ylabel('Toe Angle (degrees)')

figure('Name', 'Roll Center Height vs Roll');
plot(roll_data_ordered, rollheight_data, 'r');
title('Roll Center Height vs Roll')
xlabel('Roll Angle (degrees)')
ylabel('Roll Center Height (mm)')
