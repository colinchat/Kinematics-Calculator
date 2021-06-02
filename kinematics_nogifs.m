clear, clc, clf %, close all 

%{
    +z (up)
    |
    |
    |_______ +y (towards right)
   / 
  / 
 +x (towards back)

 * vehicle ground plane at z = 0
 * input RIGHT side of suspension (will affect toe/camber convention)
 * toe follows right hand rule about z
 * camber follows NEGATIVE right hand rule about x
 * roll coincides with wheel travel (travel up -> positive roll)
 * doesn't support inboard points changing currently (easy to add)

BONUS FUTURE FEATURES:
- lateral roll center position
- percent ackermann for steering geometry
- steering wheel to tire sensitivity 
%}

% CHANGE THESE VALUES +++++++++++++++++++++++++++++++++++++++++++++++++++++
% 21_R05 hardpoints

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
init_wheelbase = 1566.8;
CG_HEIGHT = 346; % For anti-features, from Z:\2022_master\Analysis\Suspension\Springs, Damping & Contact Patch Loads
pfb = 0.7; % percent of braking on front, rough estimate

% iteration constants
loops            = 100;    % must be even (half upwards, half downwards)
time_interval    = 0.5;    % step size in seconds 
input_vel        = 1;      % vertical wc velocity in mm/s
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% compiling starting points and saving original values
wc3_orig    = [wc_orig, wc_orig, wc_orig];
ob3_orig    = [upper, lower, tie_ob];
ob5_orig    = [ob3_orig(:,1), ob3_orig(:,1), ob3_orig(:,2), ob3_orig(:,2), ob3_orig(:,3)];
ib5_orig    = [upper_f, upper_r, lower_f, lower_r, tie_ib];
ctp_orig    = [1 0 0; ...
               0 cos(init_camber_deg*pi/180) -sin(init_camber_deg*pi/180);...
               0 sin(init_camber_deg*pi/180) cos(init_camber_deg*pi/180)]...
               * [0; 0; -1 * wc_orig(3)];

% changing variables for main loop
wc          = wc_orig;
wc3         = wc3_orig;
ob3         = ob3_orig;
ob5         = ob5_orig;
wc_to_cp   = ctp_orig;

% set up data collection arrays
vertical_data   = zeros(1, loops+1);
wheelbase_data  = zeros(1, loops+1);
htrack_data     = zeros(1, loops+1); % half of track (single wheel)
camber_data     = zeros(1, loops+1);
spin_data       = zeros(1, loops+1);
toe_data        = zeros(1, loops+1);
roll_data       = zeros(1, loops+1);
rcheight_data   = zeros(1, loops+1);
wc_to_cp_data   = zeros(3, loops+1);
ca_data_right   = zeros(15, loops+1);
ob_to_wc_data   = zeros(3, loops+1);
cp_to_ic_right  = zeros(3, loops+1);
wc_to_ic_right  = zeros(3, loops+1);
antidive_data   = zeros(1, loops+1);
antisquat_data  = zeros(1, loops+1);
% initial data points
htrack_data(1)  = init_track_mm;
camber_data(1)  = init_camber_deg;
toe_data(1)     = init_toe_deg;

% set up temporary variables to compute change
vertical_mm     = 0;
wheelbase_mm    = init_wheelbase;
track_mm        = init_track_mm; 
camber_deg      = init_camber_deg;
spin_deg        = 0;
toe_deg         = init_toe_deg;
roll_deg        = 0;

% Checks lengths of control arms at start of iteration vs end
ca5 = ib5_orig - ob5_orig;
error = [norm(ca5(:,1)), norm(ca5(:,2)), norm(ca5(:,3)), norm(ca5(:,4)), norm(ca5(:,5))];

 for iter = 1:loops
    %    Vx         Vy     Wx      Wy    Wz
    syms dWheelbase dTrack dCamber dSpin dToe
    
    % linkage vectors: 1=FUCA, 2=RUCA, 3=FLCA, 4=RLCA, 5=TROD
    ca5 = ib5_orig - ob5;
    
    % radius from wheelcenter (x3)
    rad3 = ob3 - wc3;
    
    % velocity vector at wheelcenter (x3)
    wcv3 = [dWheelbase, dWheelbase, dWheelbase; 
            dTrack,     dTrack,     dTrack;...
            input_vel,  input_vel,  input_vel];
    
    % relative velocities at outboard points
    wc_av = [dCamber; dSpin; dToe];
    wcav3 = [wc_av, wc_av, wc_av];
    obrv3 = cross(wcav3, rad3);
    
    % velocities at outboard points: 1=Upper, 2=Lower, 3=Tie
    obv3 = obrv3 + wcv3;
    
    % dot product of velocities with each control arm
    eqn1 = dot(obv3(:,1), ca5(:,1)) == 0;
    eqn2 = dot(obv3(:,1), ca5(:,2)) == 0;
    eqn3 = dot(obv3(:,2), ca5(:,3)) == 0;
    eqn4 = dot(obv3(:,2), ca5(:,4)) == 0;
    eqn5 = dot(obv3(:,3), ca5(:,5)) == 0;
    
    % solve linear equation
    sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5],...
          [dWheelbase, dTrack, dCamber, dSpin, dToe]);
    
    % subsitute variables
    obv3 = subs(obv3,...
                   [dWheelbase,     dTrack,     dCamber,     dSpin,     dToe], ...
                   [sol.dWheelbase, sol.dTrack, sol.dCamber, sol.dSpin, sol.dToe]);
    wcv = subs(wcv3(:,1),...
                   [dWheelbase,     dTrack],...
                   [sol.dWheelbase, sol.dTrack]);
    
    % integrate wheelcenter velocities
    dt = time_interval;
    ob_displ = obv3*dt; %mm/s --> mm
    wc_displ = wcv*dt;
    wc_av = [sol.dCamber; sol.dSpin; sol.dToe]; %rad/s --> rad
    wc_ang = wc_av*dt;
    
    % rotation matrices
    Rx = [1 0 0; 0 cos(wc_ang(1)) -sin(wc_ang(1)); 0 sin(wc_ang(1)) cos(wc_ang(1))];
    Ry = [cos(wc_ang(2)) 0 sin(wc_ang(2)); 0 1 0; -sin(wc_ang(2)) 0 cos(wc_ang(2))];
    Rz = [cos(wc_ang(3)) -sin(wc_ang(3)) 0; sin(wc_ang(3)) cos(wc_ang(3)) 0; 0 0 1];
    
    % rotate upright vectors, contact patch vector
    rad3 = Rz*(Ry*(Rx*rad3));
    wc_to_cp = Rx*wc_to_cp;
    
    % find vector pointing towards IC (perpendicular to V in ZY plane)
    cp_vel = cross(wc_av, wc_to_cp) + wcv;
    syms v_dist
    sol_vdist = solve(dot([0; -1; v_dist], cp_vel) == 0, v_dist);
    cp_to_ic_right(:,iter+1) = [0; -1; sol_vdist]./norm([0; -1; sol_vdist]);
    sol_vdist = solve(dot([0; -1; v_dist], wcv) == 0, v_dist);
    wc_to_ic_right(:,iter+1) = [0; -1; sol_vdist]./norm([0; -1; sol_vdist]);
    
    % record data
    wheelbase_mm    = wheelbase_mm + wc_displ(1);
    track_mm        = track_mm     + wc_displ(2);
    vertical_mm     = vertical_mm  + wc_displ(3);
    camber_deg      = camber_deg   + -wc_ang(1)*180/pi;
    spin_deg        = spin_deg     + wc_ang(2)*180/pi;
    toe_deg         = toe_deg      + wc_ang(3)*180/pi;
    
    vertical_data(iter+1)   = vertical_mm;
    wheelbase_data(iter+1)  = wheelbase_mm;
    htrack_data(iter+1)     = track_mm;
    camber_data(iter+1)     = camber_deg;
    spin_data(iter+1)       = spin_deg;
    toe_data(iter+1)        = toe_deg;
    wc_to_cp_data(:,iter+1) = wc_to_cp;
    
    % only move wheel center and rotate radius vectors
    wc  = vpa(wc + wc_displ); % vpa( ) --> round for efficiency
    wc3 = [wc, wc, wc];
    ob3 = vpa(wc3 + rad3);
    ob5 = [ob3(:,1), ob3(:,1), ob3(:,2), ob3(:,2), ob3(:,3)];
    ca5 = ib5_orig - ob5;
    ca_data_right(1:3,iter+1)     = -ca5(:,1);
    ca_data_right(4:6,iter+1)     = -ca5(:,2);
    ca_data_right(7:9,iter+1)     = -ca5(:,3);
    ca_data_right(10:12,iter+1)   = -ca5(:,4);
    ca_data_right(13:15,iter+1)   = -ca5(:,5);
    ob_to_wc_data(:,iter+1)       = -rad3(:,1);
    
    % calculate anti-features 
    [antidive_data(iter+1), antisquat_data(iter+1)] = anti_features(wc, wc_to_cp, ca5, ob5, CG_HEIGHT, wheelbase_mm, pfb);
    
    % flip from upwards input velocity to downwards at halfway
    if iter == loops/2
        fprintf('Average error in control arms: %.3fmm\n', ((norm(ca5(:,1)) - error(1))+...
            (norm(ca5(:,2)) - error(2))+...
            (norm(ca5(:,3)) - error(3))+...
            (norm(ca5(:,4)) - error(4))+...
            (norm(ca5(:,5)) - error(5)))/5);
        input_vel = -input_vel;
        wc        = wc_orig;
        wc3       = wc3_orig;
        ob5       = ob5_orig;
        ob3       = ob3_orig;
        vertical_mm    = 0;
        wheelbase_mm   = init_wheelbase;
        track_mm       = init_track_mm;
        camber_deg     = init_camber_deg;
        spin_deg       = 0;
        toe_deg        = init_toe_deg;
        roll_deg       = 0;
        wc_to_cp      = ctp_orig;
    end
 end

fprintf('Average error in control arms: %.3fmm\n', ((norm(ca5(:,1)) - error(1))+...
            (norm(ca5(:,2)) - error(2))+...
            (norm(ca5(:,3)) - error(3))+...
            (norm(ca5(:,4)) - error(4))+...
            (norm(ca5(:,5)) - error(5)))/5);

% rearranging data to go from bottom of travel to top
antidive_data   = rearrange(antidive_data, loops);
antisquat_data  = rearrange(antisquat_data, loops);
vertical_data   = rearrange(vertical_data, loops);
htrack_data     = rearrange(htrack_data, loops);
camber_data     = rearrange(camber_data, loops);
toe_data        = rearrange(toe_data, loops);
roll_data       = rearrange(roll_data, loops);
cp_to_ic_right  = rearrange(cp_to_ic_right, loops);
cp_to_ic_right(:,loops/2+1)    = cp_to_ic_right(:,loops/2);
cp_to_ic_left                  = flip(cp_to_ic_right, 2);
cp_to_ic_left(2,:)             = -1*cp_to_ic_left(2,:);
ca_data_right    = rearrange(ca_data_right, loops);
ob_to_wc_data    = rearrange(ob_to_wc_data, loops);
wc_to_cp_data    = rearrange(wc_to_cp_data, loops);
ca5_orig = ob5_orig - ib5_orig;
ca_data_right(:,loops/2+1) = [ca5_orig(:,1); ca5_orig(:,2); ca5_orig(:,3); ca5_orig(:,4); ca5_orig(:,5)];

% variables for roll
wc_to_cp_data(:,loops/2+1) = ctp_orig;
ob_to_wc_data(:,loops/2+1) = ob_to_wc_data(:,loops/2+2);
init_ground = 0;
ib5 = ib5_orig;
ib5left = [ib5(1,:); -ib5(2,:); ib5(3,:)];

for iter = loops/2+1:loops+1
    % find initial roll center
    syms rscale lscale
    rightV = [0; htrack_data(iter); 0] + (cp_to_ic_right(:,iter) .* rscale);
    leftV  = [0; -htrack_data(iter); 0] + (cp_to_ic_left(:,iter) .* lscale);
    rceqn2 = rightV(2) - leftV(2) == 0;
    rceqn3 = rightV(3) - leftV(3) == 0;
    rc_sol = solve([rceqn2, rceqn3], [rscale, lscale]);
    rc_vec = vpa(subs(rightV, rscale, rc_sol.rscale));
    rcheight_data(iter) = rc_vec(3);
    rightVisual = vpa(cp_to_ic_right(:,iter) .* rc_sol.rscale);
    leftVisual  = vpa(cp_to_ic_left(:,iter) .* rc_sol.lscale);
    rc_vec5 = [rc_vec, rc_vec, rc_vec, rc_vec, rc_vec];
    
    % rotate ob5 around roll center (increment system then make contact patch touch ground)
    syms theta
    Rx_ = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    rc_to_ib = ib5 - rc_vec5;
    rc_to_cp = rc_to_ib(:,1) + ca_data_right(1:3,iter) + ob_to_wc_data(:,iter) + wc_to_cp_data(:,iter);
    rc_to_cp = Rx_*rc_to_cp;
    cp = rc_vec + rc_to_cp;
    new_sol = 0;
    if iter == loops/2+1
        cp = subs(cp, theta, 0);
        rc_to_cp = subs(rc_to_cp, theta, 0);
        Rx_ = subs(Rx_, theta, 0);
        init_ground = cp(3);
        roll_data(iter) = 0;
    else
        rotate_eqn = cp(3) == init_ground;
        rotate_sol = solve(rotate_eqn, theta);
        [new_sol, ind] = min([abs(rotate_sol(1)), abs(rotate_sol(2))]);
        new_sol = real(rotate_sol(ind)/abs(rotate_sol(ind))*new_sol);
        rc_to_cp = subs(rc_to_cp, theta, new_sol);
        Rx_ = subs(Rx_, theta, new_sol);
        cp = rc_vec + rc_to_cp;
        roll_data(iter) = roll_data(iter-1) + vpa(new_sol)*180/pi;
    end

    %update points and vectors (relative ones)
    ib5 = rc_vec5 + Rx_*rc_to_ib; % only absolute one needs rc to not rotate
    rc_to_ib_left = ib5left - rc_vec5;
    ib5left = rc_vec5 + Rx_*rc_to_ib_left;
    ca_data_right(1:3, iter:loops+1)   = Rx_*ca_data_right(1:3, iter:loops+1);
    ca_data_right(4:6, iter:loops+1)   = Rx_*ca_data_right(4:6, iter:loops+1);
    ca_data_right(7:9, iter:loops+1)   = Rx_*ca_data_right(7:9, iter:loops+1);
    ca_data_right(10:12, iter:loops+1) = Rx_*ca_data_right(10:12, iter:loops+1);
    ca_data_right(13:15, iter:loops+1) = Rx_*ca_data_right(13:15, iter:loops+1);
    wc_to_cp_data(:,iter:loops+1)     = Rx_*wc_to_cp_data(:,iter:loops+1);
    ob_to_wc_data(:,iter:loops+1)     = Rx_*ob_to_wc_data(:,iter:loops+1);
    cp_to_ic_right(:,iter:loops+1)   = Rx_*cp_to_ic_right(:,iter:loops+1);
    cp_to_ic_left(:,iter:loops+1)    = Rx_*cp_to_ic_left(:,iter:loops+1); 
end

roll_data(1:loops/2) = -flip(roll_data(loops/2+2:loops+1), 2); % copy negative of second half
rcheight_data(1:loops/2) = flip(rcheight_data(loops/2+2:loops+1), 2); % copy second half 
rcheight_data(loops/2+1) = rcheight_data(loops/2); % fill in initial value

figure('Name', 'Camber vs Wheel Travel');
plot(vertical_data, camber_data, 'r');
title('Camber vs Wheel Travel')
xlabel('Wheel Travel (mm)')
ylabel('Camber Angle (degrees)')

figure('Name', 'Toe vs Wheel Travel')
plot(vertical_data, toe_data, 'r');
title('Toe vs Wheel Travel')
xlabel('Wheel Travel (mm)')
ylabel('Toe Angle (degrees)')

figure('Name', 'Camber vs Roll');
plot(roll_data, camber_data - roll_data, 'r');
title('Camber vs Roll')
xlabel('Roll Angle (degrees)')
ylabel('Camber Angle (degrees)')

figure('Name', 'Toe vs Roll');
plot(roll_data, toe_data, 'r');
title('Toe vs Roll')
xlabel('Roll Angle (degrees)')
ylabel('Toe Angle (degrees)')

figure('Name', 'Roll Center Height vs Roll');
plot(roll_data, rcheight_data, 'r');
title('Roll Center Height vs Roll')
xlabel('Roll Angle (degrees)')
ylabel('Roll Center Height (mm)')
set(gca, 'YLim', [0 60])

antisquat_data(loops/2+1) = (antisquat_data(loops/2)+antisquat_data(loops/2+2))/2;
figure('Name', 'Anti-squat vs Wheel travel (REAR ONLY)');
plot(vertical_data, antisquat_data, 'r');
title('Anti-squat vs Wheel travel')
xlabel('Wheel Travel (mm)')
ylabel('Percent Anti-Squat')

antidive_data(loops/2+1) = (antidive_data(loops/2)+antidive_data(loops/2+2))/2;
figure('Name', 'Anti-Dive vs Wheel travel (FRONT ONLY)');
plot(vertical_data, antidive_data, 'r');
title('Anti-dive vs Wheel travel')
xlabel('Wheel Travel (mm)')
ylabel('Percent Anti-Dive')

function x = rearrange(y, loops)
    x = [flip(y(:,loops/2+2:loops+1), 2), y(:,1:loops/2+1)];
end

function [as, ad] = anti_features(wc, wc_to_ctp, ca5, ob5, cg_height, wheelbase, pfb)
    % vertical plane through wc
    syms x y z
    pwc = y == wc(2);
    
    % upper plane
    vu1 = ca5(:,1);
    vu2 = ca5(:,2);
    posup = ob5(:,1);
    nup = cross(vu1, vu2);
    pup = dot(nup, [x;y;z]) - dot(nup, posup) == 0;
    % lower plane
    vl1 = ca5(:,3);
    vl2 = ca5(:,4);
    poslo = ob5(:,3);
    nlo = cross(vl1, vl2);
    plo = dot(nlo, [x;y;z]) - dot(nlo, poslo) == 0;
    
    % upper line
    lu = solve([pwc, pup]);
    % lower line
    ll = solve([pwc, plo]);
    
    % intersect of two lines
    ic_z = solve(ll.x == lu.x);
    ic_x = subs(ll.x, z, ic_z);
    
    % vector from contact patch
    cp = wc + wc_to_ctp;
    cp_to_svsa = [ic_x, 0, ic_z] - [cp(1), 0, cp(3)];
    theta_cp = vpa(atan(cp_to_svsa(3) / abs(cp_to_svsa(1))));
    theta_wc = vpa(atan((cp_to_svsa(3) - wc(3)) / abs(cp_to_svsa(1))));
    
    % percent anti-dive (FRONT ONLY)
    ad = pfb*tan(theta_cp)/(cg_height/wheelbase)*100;
    
    % percent anti-squat (REAR ONLY)
    as = tan(theta_wc)/(cg_height/wheelbase)*100;
end
