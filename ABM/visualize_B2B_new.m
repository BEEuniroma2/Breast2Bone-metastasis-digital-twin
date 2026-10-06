%% this main visualizes video generated



%% 
load(name);
name_video = "video_simulation\Video_B2B_"+num2str(k,"%03d")+".avi";
if not(exist("video_simulation"))
    mkdir("video_simulation");
end

w = VideoWriter(name_video,"Uncompressed AVI");
w.FrameRate = 10;
open(w)
%%
for t = 1 : length(p_sp(1).t)
    imshow(I0(:,:));

    for n_cell = 1 : length(p_sp)
        hold on

        viscircles([p_sp(n_cell).y(t),p_sp(n_cell).x(t)],20,'Color','g');
    end
    for n_cell = 1 : length(pos_ix_hbmscs)
        hold on
        viscircles([pos_iy_hbmscs(n_cell),pos_ix_hbmscs(n_cell)],30);
    end
    title("Time: "+num2str(t))
    
    
    f = getframe();
    image = imresize(f.cdata,0.2,"bilinear");
    writeVideo(w,image)
    
end

close(w)