%get_input_coords.m
%extracts relevant motion data (923,061 particles)
%from straight_guide_final_coordinates.txt

function [r_dif, v_r, y, v_y] = get_input_coords()
    fileid = fopen('straight_guide_final_coordinates.txt', 'r');
    for i=1:1:3
        fgetl(fileid);
    end

    A = textscan(fileid, '%f %f %f %f %f %f %f %d', 'Delimiter', ',');
    fclose(fileid);

    z = A{6};

    r_dif = A{2}(z >= 7.62e-4);
    v_r = A{3}(z >= 7.62e-4);
    y = A{4}(z >= 7.62e-4);
    v_y = A{5}(z >= 7.62e-4);
end
