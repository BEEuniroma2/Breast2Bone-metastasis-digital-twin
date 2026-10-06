
% file genera video B2B
clear; close all; clc;
flag_fig=false;
L = 10;
for k = 1 : L
    disp([k,L])
    name="data_simulation\video_B2B_"+num2str(k,"%03d")+".mat";
    if not(exist("data_simulation"))
        mkdir("data_simulation");
    end
    
    settings_file="settings_file_B2B.m";
    tic
    my_cell_motility_model(name,settings_file,flag_fig);
    
    toc
    run("visualize_B2B_new.m")
end


