function coop_11_4()
    %%% change these!
    ti = 0;
    tf = 10;
    %struct values
    m = 1;
    I = 1;
    g = -9.81;
    ks = [10;10;10;10;10;10;10;10;];
    l0s = [0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;];
    Pw = [1,1;-1,1;-1.5,0.5;-1.5,-0.5;-1,-1;1,-1;1.5,-0.5;1.5,0.5;]';
    Pb = [1,0.5;-1,0.5;-1,0.5;-1,-0.5;-1,-0.5;1,-0.5;1,-0.5;1,0.5;]';
    %start condition
    x0 = 0;
    y0 = -0.2253;
    theta0 = 0;
    vx0 = 0;
    vy0 = 0;
    vtheta0 = 0;

    %%%don't change these!
    %set up for integration
    box_params = struct_maker(m,I,g,ks,l0s,Pw,Pb);
    my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
    Veq_approx = [x0;y0;theta0;vx0;vy0;vtheta0];
    tspan = [ti,tf];
    DormandPrince = makedp();

    %linearized sim
    eqFUN = @(V_in) box_rate_func(0,V_in,box_params);
    [Veq, it, flag, glist] = multi_newton_solver(eqFUN, Veq_approx, 1e-14, 1e-14, 200, 1);
    J_approx = approximate_jacobian(eqFUN, Veq);
    Q = - J_approx(4:6,1:3);
    my_linear_rate = @(t_in,V_in) J_approx*(V_in-V_eq);

    epsilon = 1e-1;
    V0 = Veq + epsilon*Veq_approx;
    
    %integration
    h0 = 0.1;
    des_err = 10^-2.5;
    [t_list,X_list,h_avg, num_evals, fail_rate, h_rec] = explicit_RK_variable_step_integration(my_rate_func, tspan, V0, h0, DormandPrince, 5, des_err);
    plot(t_list, X_list(2,:))

    [V_mat,D_mat] = eig(Q);

    mode_index = 3;
    %mode shape
    U_mode = V_mat(:,mode_index);
    omega_n = sqrt(-D_mat(mode_index,mode_index));
    %small number
    epsilon = 10e-14;
    V0 = Veq + epsilon*[U_mode;0;0;0];
    %run the integration of nonlinear system
    [t_list,X_list,h_avg, num_evals, fail_rate, h_rec] = explicit_RK_variable_step_integration(my_rate_func, tspan, V0, h0, DormandPrince, 5, des_err);

    x_modal = Veq(1)+epsilon*U_mode(1)*cos(omega_n*t_list);
    y_modal = Veq(2)+epsilon*U_mode(2)*cos(omega_n*t_list);
    theta_modal = Veq(3)+epsilon*U_mode(3)*cos(omega_n*t_list);
    figure
    hold on
    plot(t_list,y_modal)
    plot(t_list, X_list(2,:))
end

function box_params = struct_maker(m,I,g,ks,l0s,Pw,Pb)
    box_params = struct();
    box_params.m = m;
    box_params.I = I;
    box_params.g = g;
    box_params.k_list = ks;
    box_params.l0_list = l0s;
    box_params.P_world = Pw;
    box_params.P_box = Pb;
end

function F = compute_spring_force(k,l0,PA,PB)
    %current length of the spring
    l = norm(PB-PA);
    %unit vector pointing from PA to PB
    e_s = (PB-PA)/l;
    %Force exerted by spring at point B
    F = -k.*(l-l0).*e_s;
end

function Plist_world = compute_rbt(x,y,theta,Plist_box)
    rc = [x;y;];
    R = [cosd(theta), -sind(theta); sind(theta), cosd(theta);];
    [~,N] = size(Plist_box);
    offset = repmat(rc,1,N);
    Plist_world = R*Plist_box+offset;
end

function [ax,ay,atheta] = compute_accel(x,y,theta,box_params)
    rc = [x;y;0];
    m = box_params.m;
    ks = box_params.k_list;
    g = box_params.g;
    l0s = box_params.l0_list;
    num_s = length(ks);
    Fs = zeros(2, num_s);
    Ts = zeros(1, num_s);
    Pb_world = compute_rbt(x,y,theta,box_params.P_box);
    Pm_world = box_params.P_world;
    for i = 1:num_s
        Pm = [Pm_world(:,i);0];
        Pb = [Pb_world(:,i);0];
        F = compute_spring_force(ks(i),l0s(i),Pm,Pb);
        T = cross((Pb-rc), F);
        Fs(:, i) = F(1:2);
        Ts(i) = T(3);
    end
    gvec = [0;g];
    axy = gvec+1/m.*sum(Fs,2);
    ax = axy(1);
    ay = axy(2);
    atheta = sum(Ts)/box_params.I;
end

function dVdt = box_rate_func(t,V,box_params)
    x=V(1); y=V(2); theta=V(3);
    dxdt=V(4); dydt=V(5); dthetadt=V(6);
    [ax,ay,atheta] = compute_accel(x,y,theta,box_params);
    dVdt = [dxdt;dydt;dthetadt;ax;ay;atheta;];
end


%%%numerical integrator
function [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct)
    as = BT_struct.A;
    bs = BT_struct.B;
    cs = BT_struct.C;
    s = length(cs);
    m = length(XA);
    ks = zeros(m, s);
    for i = 1:s
        SUMak = ks*(as(i,:)');
        ki = rate_func_in((t+cs(i)*h), (XA+h*SUMak));
        ks(:, i) = ki;
    end
    SUMbk1 = ks*bs(1,:)';
    SUMbk2 = ks*bs(2,:)';
    XB1 = XA + h*SUMbk1;
    XB2 = XA + h*SUMbk2;
    num_evals = s;
end

function [XB, num_evals, h_next, redo] = explicit_RK_variable_step(rate_func_in,t,XA,h,BT_struct,p,error_desired)
    alpha = 1.5;

    redo = false;
    [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct);
    epsc = norm(XB1-XB2);
    temp = (error_desired/epsc).^(1/p);
    h_next = min(0.9*temp, alpha)*h;
    XB = XB1;
    if error_desired < epsc
        redo = true;
    end
end

function [t_list, X_list, h_avg, num_fails, num_evals, h_rec] = variable_step_integration(rate_func_in, step_func, tspan, X0, h_ref, BT_struct, p, error_desired)
    ti = tspan(1); tf = tspan(2);
    N  = ceil((tf - ti)/h_ref); h = (tf - ti)/N;
    nx = numel(X0); 
    t_list = 1; t_list(1) = ti;
    X_list = zeros(nx, 1); X_list(:,1) = X0;
    num_evals = 0; XA = X0;
    h_rec = []; 
    num_fails = 0;
    tnow = ti;
    while tnow <= tf
        redo = true; 
        while redo == true
            h_prev = h;
            h = min(h,tf-tnow+1e-15);
            [XB, add_evals, h, redo] = step_func(rate_func_in, tnow, XA, h, BT_struct, p, error_desired);
            if redo == true
                num_fails = num_fails+1;
            end
            num_evals = num_evals + add_evals;
        end
        h_rec(end+1,:) = h_prev;
        tnow = tnow+h_prev;
        t_list(end+1) = tnow;
        X_list(:,end+1) = XB; XA = XB;
    end
    h_avg = mean(h_rec);
end

function [t_list,X_list,h_avg, num_evals, fail_rate, h_rec] = explicit_RK_variable_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct, p, error_desired)
    [t_list,X_list,h_avg, num_fails, num_evals, h_rec] = variable_step_integration(rate_func_in, @explicit_RK_variable_step, tspan, X0, h_ref, BT_struct, p, error_desired);
    fail_rate = num_fails/num_evals;
end


%%%imported functions
% jacobian
function J = approximate_jacobian(FUN, X)
% Centered finite-difference Jacobian: J(:,i) ≈ (F(x+h/2*ei)-F(x-h/2*ei))/h
    h = 1e-6;
    n = numel(X);
    F0 = FUN(X);
    m = numel(F0);
    J = zeros(m, n);
    for i = 1:n
        ei      = zeros(n,1); ei(i) = h/2;
        F_plus  = FUN(X + ei);
        F_minus = FUN(X - ei);
        J(:,i)  = (F_plus - F_minus) / h;
    end
end

%multi newton solver
function [root, it, flag, glist] = multi_newton_solver(FUN_both, X0, Athresh, Bthresh, maxit, numdiff)
% General multidimensional Newton’s method.
% FUN_both(X): either returns F(X)            (numdiff==1)   OR
%                 returns [F(X), J(X)]          (numdiff==2)
% Early termination: ||ΔX|| < Athresh OR ||F|| < Bthresh
% Safeguard: stop if Jacobian is (near) singular
% Outputs:
%     root: solution vector
%     it  : iteration count
%     flag: 1 success, -2 bad/singular Jacobian
%     glist: columns are the iterate history
    glist = [];
    root  = X0;
    it    = 0;
    flag  = 0;
    glist(:,end+1) = root;
    % get residual & Jacobian
    if numdiff == 2
        try
            [FX, J] = FUN_both(root);
        catch
            FX = FUN_both(root);
            J  = approximate_jacobian(FUN_both, root);
        end
    else
        FX = FUN_both(root);
        J  = approximate_jacobian(FUN_both, root);
    end
    while norm(FX) > Bthresh && it < maxit
        % singularity guard (use J*J' to handle non-square cases robustly)
        if abs(det(J*J')) < eps
            flag = -2; return;
        end
        X_new = root - J\FX;
        glist(:,end+1) = X_new;
        if norm(X_new - root) < Athresh
            root = X_new; flag = 1; return
        end
        root = X_new;
        if numdiff == 2
            try
                [FX, J] = FUN_both(root);
            catch
                FX = FUN_both(root);
                J  = approximate_jacobian(FUN_both, root);
            end
        else
            FX = FUN_both(root);
            J  = approximate_jacobian(FUN_both, root);
        end
        it = it + 1;
    end
    flag = 1; % reached residual threshold or maxit
end

%%%get dormand prince out of the way
function DormandPrince = makedp()
    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
    1/5, 0, 0, 0,0,0,0;...
    3/40, 9/40, 0, 0, 0, 0,0;...
    44/45, -56/15, 32/9, 0, 0, 0,0;...
    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];
end


%%%ANIMATION
%====================================================================
% AVI recorder: simple box + straight springs animation
%====================================================================
function record_animation_avi(filename, frame_rate, t_list, X_list, box_params)

    v = VideoWriter(filename, 'Motion JPEG AVI');
    v.FrameRate = frame_rate;
    open(v);

    fig = figure('Color','white');
    axis equal; axis([-3 3 -3 3]); hold on;
    xlabel('x'); ylabel('y');
    title('vibrating box animation');

    % square box corners in box frame (side length 2)
    box_corners_box = [ 1  -1  -1   1   1;  ...
                        1   1  -1  -1   1];

    Pm_world = box_params.P_world;   % fixed anchors (2 x n)

    for k = 1:length(t_list)
        x     = X_list(1,k);
        y     = X_list(2,k);
        theta = X_list(3,k);

        % box corners and spring attachment points in world frame
        Pb_world = compute_rbt(x,y,theta, box_params.P_box);
        box_corners_world = compute_rbt(x,y,theta, box_corners_box);

        cla; hold on; axis equal; axis([-3 3 -3 3]);

        % draw outer frame (just a big square for context)
        rectangle('Position',[-2.5 -2.5 5 5], ...
                  'EdgeColor',[0.7 0.7 0.7]);

        % draw springs as straight lines
        for i = 1:size(Pb_world,2)
            plot([Pm_world(1,i) Pb_world(1,i)], ...
                 [Pm_world(2,i) Pb_world(2,i)], 'k-');
            plot(Pm_world(1,i), Pm_world(2,i), 'ro', ...
                 'MarkerFaceColor','r','MarkerSize',4);
            plot(Pb_world(1,i), Pb_world(2,i), 'ro', ...
                 'MarkerFaceColor','r','MarkerSize',4);
        end

        % draw the box
        patch(box_corners_world(1,:), box_corners_world(2,:), ...
              [0.8 0.8 0.8], 'EdgeColor','k');

        title(sprintf('vibrating box, t = %.2f s', t_list(k)));

        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end

    close(v);
    close(fig);

    fprintf('AVI animation saved to %s\n', filename);
end

%====================================================================
% spring plotting utilities (from handout)
%====================================================================
function spring_plot_struct = initialize_spring_plot(num_zigs,w)
    spring_plot_struct = struct();
    zig_ending = [.25,.75,1; ...
                  -1,1,0];
    zig_zag = zeros(2,3+3*num_zigs);
    zig_zag(:,1)   = [-.5;0];
    zig_zag(:,end) = [num_zigs+.5;0];
    for n = 0:(num_zigs-1)
        zig_zag(:,(3+3*n):2+3*(n+1)) = zig_ending + [n,n,n;0,0,0];
    end
    zig_zag(1,:) = (zig_zag(1,:)-zig_zag(1,1))/(zig_zag(1,end)-zig_zag(1,1));
    zig_zag(2,:) = zig_zag(2,:)*w;
    spring_plot_struct.zig_zag    = zig_zag;
    spring_plot_struct.line_plot  = plot(0,0,'k','linewidth',2);
    spring_plot_struct.point_plot = plot(0,0,'ro','markerfacecolor','r','markersize',7);
end

function update_spring_plot(spring_plot_struct,P1,P2)
    dP = P2-P1;
    R = [dP(1), -dP(2)/norm(dP); ...
         dP(2),  dP(1)/norm(dP)];
    plot_pts = R*spring_plot_struct.zig_zag;
    set(spring_plot_struct.line_plot,...
        'xdata',plot_pts(1,:)+P1(1),...
        'ydata',plot_pts(2,:)+P1(2));
    set(spring_plot_struct.point_plot,...
        'xdata',[P1(1),P2(1)],...
        'ydata',[P1(2),P2(2)]);
end