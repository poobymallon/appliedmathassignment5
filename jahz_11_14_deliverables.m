function jahz_11_14_deliverables()
% main script for assignment 5 deliverables
% generates:
%  1) centroid path and printed equilibrium
%  2) linear vs nonlinear comparison (small and large perturbations)
%  3) modal analysis plots for 3 modes
%  4) an .avi animation of the box vibrating in mode 1

    %-----------------------------
    % system parameters and RK data
    %-----------------------------
    box = make_box_params();        % struct with geometry and springs
    DP  = make_DP_tableau();        % Dormand–Prince 5(4) embedded RK
    p_ord   = 5;                    % order of higher method in pair
    h0      = 0.05;                 % initial step guess
    err_des = 1.0e-5;               % desired local error
    tspan   = [0 10];               % 10 seconds of motion

    %-----------------------------
    % 1) find an equilibrium
    %-----------------------------
    V_guess = [0; -0.2; 0; 0; 0; 0];     % [x;y;theta;vx;vy;vtheta]
    [V_eq, newton_hist] = newton_equilibrium(@(V) box_rate_func(0,V,box), V_guess);

    fprintf('equilibrium found at [x y th] = [%.4f %.4f %.4f]\n', ...
        V_eq(1), V_eq(2), V_eq(3));

    % quick centroid path from some offset around equilibrium
    V0_centroid = V_eq + [0.1; 0.05; 0.15; 0; 0; 0];
    [t_centroid, V_centroid] = rk_adaptive(@(t,V) box_rate_func(t,V,box), ...
        tspan, V0_centroid, h0, DP, p_ord, err_des);

    figure; hold on; grid on; axis equal;
    plot(V_centroid(1,:), V_centroid(2,:), 'k-');
    plot(V_eq(1), V_eq(2), 'ko', 'markerfacecolor','k');
    xlabel('x'); ylabel('y');
    title('centroid path');
    text(0.4,0.25, sprintf('equilibrium found at [x y th] = [%.4f %.4f %.4f]', ...
        V_eq(1), V_eq(2), V_eq(3)));

    %-----------------------------
    % 2) linearization at equilibrium
    %-----------------------------
    A = J_approx(@(V) box_rate_func(0,V,box), V_eq);   % 6x6 Jacobian at eq

    % linear rate: d/dt ΔV = A (V - V_eq)
    lin_rate = @(t,V) A*(V - V_eq);

    %-----------------------------
    % 3) linear vs nonlinear comparison
    %-----------------------------
    eps_small = 0.02;   % small perturbation
    eps_large = 0.2;    % large perturbation

    dir = [1;0;0;0;0;0];   % bump in x only

    % small
    V0_sm = V_eq + eps_small*dir;
    [t_nl_sm, V_nl_sm] = rk_adaptive(@(t,V) box_rate_func(t,V,box), ...
        tspan, V0_sm, h0, DP, p_ord, err_des);
    [t_li_sm, V_li_sm] = rk_adaptive(lin_rate, ...
        tspan, V0_sm, h0, DP, p_ord, err_des);

    % large
    V0_lg = V_eq + eps_large*dir;
    [t_nl_lg, V_nl_lg] = rk_adaptive(@(t,V) box_rate_func(t,V,box), ...
        tspan, V0_lg, h0, DP, p_ord, err_des);
    [t_li_lg, V_li_lg] = rk_adaptive(lin_rate, ...
        tspan, V0_lg, h0, DP, p_ord, err_des);

    % plot small (Δx)
    figure; hold on; grid on;
    dx_nl = V_nl_sm(1,:) - V_eq(1);
    dx_li = V_li_sm(1,:) - V_eq(1);
    plot(t_nl_sm, dx_nl, 'k-', 'displayname','nonlinear x');
    plot(t_li_sm, dx_li, 'r--', 'displayname','linear x');
    xlabel('t'); ylabel('\Delta x');
    title('linear vs nonlinear from small perturbation');
    legend('location','best');

    % plot large (Δx)
    figure; hold on; grid on;
    dx_nl_L = V_nl_lg(1,:) - V_eq(1);
    dx_li_L = V_li_lg(1,:) - V_eq(1);
    plot(t_nl_lg, dx_nl_L, 'k-', 'displayname','nonlinear x');
    plot(t_li_lg, dx_li_L, 'r--', 'displayname','linear x');
    xlabel('t'); ylabel('\Delta x');
    title('linear vs nonlinear from large perturbation');
    legend('location','best');

    %-----------------------------
    % 4) modal analysis
    %-----------------------------
    % A has block structure [0 I; -Q 0]
    Q = -A(4:6,1:3);
    [U_modes, Lambda] = eig(Q);
    wn = sqrt(diag(Lambda));      % resonant frequencies (rad/s)

    % sort by frequency
    [wn, idx] = sort(wn);
    U_modes = U_modes(:,idx);

    fprintf('\nresonant frequencies (rad/s)\n');
    for i = 1:3
        fprintf('  mode %d:  omega = %.3f\n', i, wn(i));
    end

    eps_mode = 0.1;
    t_mode = t_nl_sm;     % reuse same t-grid length via nonlinear sim

    figure;
    for mode_idx = 1:3
        U0 = U_modes(:,mode_idx);
        omega_n = wn(mode_idx);

        % nonlinear sim with mode-type initial displacement
        V0_mode = V_eq + eps_mode*[U0;0;0;0];
        [t_nl, V_nl] = rk_adaptive(@(t,V) box_rate_func(t,V,box), ...
            tspan, V0_mode, h0, DP, p_ord, err_des);

        % modal prediction: ΔU(t) = eps_mode * U0 * cos(omega_n t)
        t_modal = t_nl;
        cos_term = cos(omega_n * t_modal);
        dx_modal = eps_mode * U0(1) * cos_term;
        dy_modal = eps_mode * U0(2) * cos_term;
        dth_modal= eps_mode * U0(3) * cos_term;

        dx_nl = V_nl(1,:) - V_eq(1);
        dy_nl = V_nl(2,:) - V_eq(2);
        dth_nl= V_nl(3,:) - V_eq(3);

        subplot(3,1,mode_idx); hold on; grid on;
        plot(t_nl, dx_nl, 'k-', 'displayname','nonlinear \Deltax');
        plot(t_modal, dx_modal, 'r--', 'displayname','modal \Deltax');
        plot(t_nl, dy_nl, 'g-', 'displayname','nonlinear \Deltay');
        plot(t_modal, dy_modal, 'c--', 'displayname','modal \Deltay');
        plot(t_nl, dth_nl, 'b-', 'displayname','nonlinear \Delta\theta');
        plot(t_modal, dth_modal,'m--', 'displayname','modal \Delta\theta');
        xlabel('t');
        ylabel('displacement');
        title(sprintf('mode %d   \\omega = %.3f rad/s', mode_idx, omega_n));
        if mode_idx == 1
            legend('location','best');
        end
    end

    %-----------------------------
    % 5) create an .avi animation for mode 1
    %-----------------------------
    mode1 = U_modes(:,1);
    V0_anim = V_eq + eps_mode*[mode1;0;0;0];
    [t_anim, V_anim] = rk_adaptive(@(t,V) box_rate_func(t,V,box), ...
        tspan, V0_anim, h0, DP, p_ord, err_des);

    record_animation_avi('mode1_box.avi', 30, t_anim, V_anim, box);
end

%====================================================================
% helper: system parameters
%====================================================================
function box = make_box_params()
    % eight springs attaching a square box to a square frame
    m = 1;
    I = 1;
    g = -9.81;

    ks  = 10*ones(8,1);
    l0s = 0.5*ones(8,1);

    % world-frame static anchor points (columns)
    Pw = [ 1,   1,  -1.5,-1.5,-1,  1,  1.5, 1.5; ...
           1,  -1,   0.5,-0.5,-1,-1,-0.5, 0.5];

    % box-frame mounting points (columns)
    Pb = [ 1,  1, -1, -1, -1,  1,  1,  1; ...
           0.5,-0.5, 0.5,-0.5,-0.5,-0.5,-0.5,0.5];

    box = struct();
    box.m        = m;
    box.I        = I;
    box.g        = g;
    box.k_list   = ks(:);
    box.l0_list  = l0s(:);
    box.P_world  = Pw;
    box.P_box    = Pb;
end

%====================================================================
% helper: Dormand–Prince 5(4) Butcher tableau
%====================================================================
function DP = make_DP_tableau()
    DP.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DP.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0; ...
            5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DP.A = [0,0,0,0,0,0,0; ...
            1/5,0,0,0,0,0,0; ...
            3/40,9/40,0,0,0,0,0; ...
            44/45,-56/15,32/9,0,0,0,0; ...
            19372/6561,-25360/2187,64448/6561,-212/729,0,0,0; ...
            9017/3168,-355/33,46732/5247,49/176,-5103/18656,0,0; ...
            35/384,0,500/1113,125/192,-2187/6784,11/84,0];
end

%====================================================================
% spring force: 3D version (z component 0)
%====================================================================
function F = compute_spring_force(k,l0,PA,PB)
    % PA, PB 3x1 vectors
    r = PB - PA;
    l = norm(r);
    e_s = r / l;
    F = -k*(l - l0)*e_s;
end

%====================================================================
% rigid body transform from box frame to world frame for points
%====================================================================
function Plist_world = compute_rbt(x,y,theta,Plist_box)
    R = [cos(theta), -sin(theta); ...
         sin(theta),  cos(theta)];
    rc = [x; y];
    [~,N] = size(Plist_box);
    Plist_world = R*Plist_box + rc*ones(1,N);
end

%====================================================================
% acceleration at pose (x,y,theta)
%====================================================================
function [ax,ay,atheta] = accel_at_pose(x,y,theta,box)
    rc       = [x;y;0];
    m        = box.m;
    I        = box.I;
    ks       = box.k_list;
    l0s      = box.l0_list;
    num_s    = numel(ks);

    Pb_world = compute_rbt(x,y,theta, box.P_box);
    Pm_world = box.P_world;

    Fs = zeros(2,num_s);
    Ts = zeros(1,num_s);

    for i = 1:num_s
        Pm = [Pm_world(:,i);0];
        Pb = [Pb_world(:,i);0];
        F  = compute_spring_force(ks(i), l0s(i), Pm, Pb);
        tau = cross(Pb-rc, F);
        Fs(:,i) = F(1:2);
        Ts(i)   = tau(3);
    end

    gvec = [0; box.g];
    axy  = gvec + (1/m)*sum(Fs,2);
    ax   = axy(1);
    ay   = axy(2);
    atheta = sum(Ts)/I;
end

%====================================================================
% rate function for full box dynamics
%====================================================================
function dVdt = box_rate_func(~,V,box)
    x = V(1);  y = V(2);  theta = V(3);
    vx = V(4); vy = V(5); vtheta = V(6);

    [ax,ay,atheta] = accel_at_pose(x,y,theta,box);
    dVdt = [vx; vy; vtheta; ax; ay; atheta];
end

%====================================================================
% embedded RK step for adaptive DP
%====================================================================
function [XB1, XB2, num_evals] = rk_step_embedded(rate,t,XA,h,BT)
    A = BT.A;
    B = BT.B;
    C = BT.C;
    s = numel(C);
    n = numel(XA);
    K = zeros(n,s);
    for i = 1:s
        a = A(i,1:i-1);
        sum_prev = K(:,1:i-1)*a';
        K(:,i) = rate(t + C(i)*h, XA + h*sum_prev);
    end
    XB1 = XA + h*(K*B(1,:)');
    XB2 = XA + h*(K*B(2,:)');
    num_evals = s;
end

function [XB, num_evals, h_next, redo] = rk_step_adaptive(rate,t,XA,h,BT,p,err_des)
    alpha = 1.5;
    redo = false;

    [XB1, XB2, num_evals] = rk_step_embedded(rate,t,XA,h,BT);

    eps_c = norm(XB1 - XB2);
    temp  = (err_des/eps_c)^(1/p);
    h_next = min(0.9*temp, alpha)*h;
    XB = XB1;

    if err_des < eps_c
        redo = true;
    end
end

%====================================================================
% variable-step integration wrapper
%====================================================================
function [t_list,X_list,h_avg,num_fails,num_evals,h_rec] = ...
    rk_variable(rate,tspan,X0,h_ref,BT,p,err_des)

    ti = tspan(1); tf = tspan(2);
    N  = ceil((tf - ti)/h_ref);
    h  = (tf - ti)/N;

    t_list = ti;
    X_list = X0;
    h_rec  = [];
    num_fails = 0;
    num_evals = 0;

    t_now = ti;
    XA = X0;

    while t_now < tf - 1e-14
        redo = true;
        while redo
            h_prev = h;
            h = min(h, tf - t_now);
            [XB, adds, h, redo] = rk_step_adaptive(rate, t_now, XA, h, BT, p, err_des);
            if redo
                num_fails = num_fails + 1;
            end
            num_evals = num_evals + adds;
        end

        h_rec(end+1,1) = h_prev;
        t_now = t_now + h_prev;
        t_list(end+1,1) = t_now;
        X_list(:,end+1) = XB;
        XA = XB;
    end

    h_avg = mean(h_rec);
end

function [t_list,X_list] = rk_adaptive(rate,tspan,X0,h_ref,BT,p,err_des)
    [t_list,X_list,~,~,~,~] = rk_variable(rate,tspan,X0,h_ref,BT,p,err_des);
end

%====================================================================
% Newton solver for equilibrium f(V) = 0
%====================================================================
function [V_eq, hist] = newton_equilibrium(f,V0)
    max_iter = 20;
    tol_step = 1e-10;

    V = V0;
    hist = zeros(max_iter,1);
    for k = 1:max_iter
        F = f(V);
        J = J_approx(f,V);
        dV = -J\F;
        V  = V + dV;
        hist(k) = norm(F);
        if norm(dV) < tol_step
            hist = hist(1:k);
            break;
        end
    end
    V_eq = V;
end

% finite-difference Jacobian
function J = J_approx(f,V)
    n = numel(V);
    J = zeros(n);
    h = 1e-6;
    f0 = f(V);
    for j = 1:n
        e = zeros(n,1); e(j) = 1;
        fj = f(V + h*e);
        J(:,j) = (fj - f0)/h;
    end
end

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

