%get_survival_rate.m

%     parameters                     type
%r_dif, v_r, v_y, v_long            gpuArray
%v_long, R_tube, r_tube, m, dt      double

% output: survival_rate             double
% percent of particles that don't hit the boundary by pi/2 rad along path

function survival_rate = get_survival_rate(r_dif, v_r, y, v_y, v_long, R_tube, r_tube, m, dt)
    n_total_particles = numel(r_dif);
    len = ceil(0.5*pi*(R_tube + r_tube)/v_long/dt);
    coefs = [-2.13226e-3, 2.18636e-6, 5.92449e-8, -2.16494e-10, 2.56143e-13, 2.89048e-16];
    coder.gpu.constantMemory(coefs);

    r = r_dif + R_tube;
    theta = gpuArray.zeros(n_total_particles, 1);
    omega = v_long./r;
    l = r.*v_long;

    d = max(hypot(r - R_tube, y), eps);
    f = -sign(d).*polyval(coefs, abs(d)).*abs(d).^2;
    a_r = f.*(r - R_tube)./d/m + omega.^2.*r;
    a_y = f.*y./d/m;
    
    alv = gpuArray(true(n_total_particles, 1));
    num_dead = 0;

    dt2 = dt^2;
    half_dt = 0.5*dt;
    half_pi = 0.5*pi;
    
    for iteration = 1:len
        r = r + (v_r*dt + 0.5*a_r*dt2).*alv;
        y = y + (v_y*dt + 0.5*a_y*dt2).*alv;
        d = max(hypot(r - R_tube, y).*alv, eps);

        omega_old = omega;
        omega = l./r.^2.*alv;
        theta = theta + half_dt.*(omega_old + omega).*alv;

        new_dead = alv & (d >= r_tube);
        finished = alv & (theta >= half_pi);
        alv(new_dead | finished) = false;

        num_dead = num_dead + sum(new_dead);
        if nnz(alv) == 0
            break;
        end

        a_r_old = a_r;
        a_y_old = a_y;

        f = -sign(d).*polyval(coefs, abs(d)).*d.^2.*alv;
        a_r = f.*(r - R_tube)./d/m + omega.^2.*r;
        a_y = f.*y./d/m;

        v_r = v_r + half_dt*(a_r_old + a_r).*alv;
        v_y = v_y + half_dt*(a_y_old + a_y).*alv;
    end

    survival_rate = (1 - num_dead/n_total_particles)*100;
end
