function [t_values, b_values] = deterministicPointModel(N_i, input_time_points, STOP_TIME, dt, distance, D, k_on, k_off)

    t_values = (0:dt:STOP_TIME).';
    b_values = zeros(size(t_values));

    for idx = 2:length(t_values)

        t_now = t_values(idx);

        c_t = 0;

        for pulseIdx = 1:length(N_i)

            tau = t_now - input_time_points(pulseIdx);

            if tau > 0
                h = 1 ./ ((4*pi*D*tau).^(3/2)) .* exp(-distance^2 ./ (4*D*tau));
                c_t = c_t + N_i(pulseIdx) * h;
            end
        end

        b_prev = b_values(idx-1);
        db_dt = k_on * c_t * (1 - b_prev) - k_off * b_prev;

        b_values(idx) = b_prev + dt * db_dt;
        b_values(idx) = max(0, min(1, b_values(idx)));
    end
end