function my_cell_motility_model(name_save,settings_file,flag_fig)
% This function implements the model for Breast to Bone metastasis on chip
% simulation  
% INPUT Arguments
% - name_save: name of the file to be generated where the data will be
%              stored
% - settings_file: file where parameters of the simulation are taken
% - flag_fig:  flag for visulazation
%   Detailed explanation goes here
    visualize = false;
    %% first get parameters from the settings file
    run(settings_file);
    load(settings_file{1}(1:end-1)+"mat",'beta','a0','a1','a0_hbmscs','a1_hbmscs','alpha','as','cell_dim','D0','D_lead','n_cancer','n_hbmscs','perc_lead','R_int_c','R_int_c_far','r_lead','r_surface','r_c','r_hbmscs','rho0','sigma_lead','sigmax_0','sigmax_1','sigmay_0','sigmay_1','t0','U0','U0_far','a0_far','U1','Ul','U0_hbmscs','U1_hbmscs','Us','vmax','D0_x','D0_y');
    
    % fix the number of frames 
    N_frames = 100;
    % starting time for the cancer cells (the hbmscs are fixed in the space)
    Tstart = ones(n_cancer,1);
    
    % Get the spatial dimension of the image
    Nx = 1024;
    Ny = 1024;
    
    % generate the grid where cells are able to move
    
    x = 1 : Nx;
    y = 1 : Ny;
    [Y,X] = meshgrid(y,x);
    % initialize positions of cells
    I0 = 0*ones(Nx,Ny,1,'uint8');
    flag(1:n_cancer) = false;
    % pos_ix = r_c + randi([1 round(size(I0,1))],n_cancer,1)+Nx; % la riga deve spaziare sull'intera colonna
    pos_ix = r_c + randi([1 round(size(I0,1))],n_cancer,1);
    % 
    % pos_ix=linspace(r_sp,Nx-r_sp,n_sp);
    % pos_iy = -2*r_c + (3*r_c+2*r_c)*rand(n_cancer,1)+Ny;
    pos_iy = randn(n_cancer,1)*Ny/5;
    % pos_iy = randi([-round(4*r_c), round(4*r_c)],n_cancer,1);

    % initialize position for human bone marrow stromal cells
    % pos_ix_hbmscs = size(I0,1)/2;
    % pos_iy_hbmscs = size(I0,2)/2;

    % pos_ix = size(I0,1)/2;
    % pos_iy = size(I0,2)/4;
    pos_ix_hbmscs = randi([round(r_hbmscs), round(Nx)-round(r_hbmscs)],n_hbmscs,1);

    pos_iy_hbmscs = randi([round(r_hbmscs), round(Ny)-round(r_hbmscs)],n_hbmscs,1);

    % hbmscs are not going to change throughout the whole video
    U0_ = rand(n_cancer,1)*U0;
    U1_ = rand(n_cancer,1)*U1;
    U0_hbmscs_ = (1+(randn(n_cancer,1)/3))*U0_hbmscs;
    U1_hbmscs_ = (1+(randn(n_cancer,1)/3))*U1_hbmscs;
    


    % U0_hbmscs_ = rand(n_cancer,1)*U0_hbmscs;
    % U1_hbmscs_ = rand(n_cancer,1)*U1_hbmscs;
    % U0_far_ = (rand(n_cancer,1)-0.5)*U0_far;
    v1 = randn(n_cancer/2,1)*U0_far/10;
    v2 = (1+randn(n_cancer/2,1))*U0_far;
    v = [v1;v2];
    U0_far_ = v(randperm(length(v)));
    theta_far_ = randn(n_cancer,1)*pi/6;%
    % U0_far_ = (randn(n_cancer,1))*U0_far;
    

    % get initial positions in cancer (moving) cells;
    for k = 1 : n_cancer
        k;
        p_sp_start(k).y=pos_iy(k);
        p_sp_start(k).x=pos_ix(k);
        idx2(k).vec=[];
        idx3(k).vec=[];
    end

    for k = 1 : n_cancer
        flag(k) = true;
        p_sp(k).x(1) =  p_sp_start(k).x;
        p_sp(k).y(1) =  p_sp_start(k).y;
        p_sp(k).t(1) = 1;
        vel_x(k,1)=0;
        vel_y(k,1)=0;
        dx_nonint(k) = 0; %deviazione media del primo passo in x (asse verticale)
        dy_nonint(k) = D0; %deviazione media del primo passo in y (asse orizzontale)
        p_sp(k).x(2) =  p_sp(k).x(1) + dx_nonint(k)+sigmax_0*randn();
        p_sp(k).y(2) =  p_sp(k).y(1) + dy_nonint(k)+sigmay_0*randn();
        p_sp(k).t(2) = 2;
        % questa è la versione NCA
        % vel_x(k,2) = 0;
        % vel_y(k,2) = 0;
        % questa è la mia proposta
        vel_x(k,2) = p_sp(k).x(2)-p_sp(k).x(1);

        vel_y(k,2) = p_sp(k).y(2)-p_sp(k).y(1);
    end

    D0_x_ = randn(n_cancer,1)*D0_x;
    D0_y_ = randn(n_cancer,1)*D0_y;
    % start the loop over time
    for i = 3 : N_frames
        if mod(i,12)==0
            U0_far_ = U0_far_(randperm(length(U0_far_)));
            D0_x_ = D0_x_(randperm(length(D0_x_)));
            D0_y_ = D0_y_(randperm(length(D0_y_)));


        end
        % disp([i, N_frames])
        for k = 1 : n_cancer
            p_cx(k,i-1) = min(3*Nx,p_sp(k).x(i-1));
            p_cy(k,i-1) = min(3*Ny,p_sp(k).y(i-1));
        end


        %% loop over cells for analyzing the 
        for k = 1 : n_cancer

            % initialize variables

            % % verifica se è leader
            % xn=0;
            % yn=0;
            % xn_unif_v=0;
            % yn_unif_v=0;
            % xs=0;
            % ys=0;
            % xl=0;
            % yl=0;
            % xs2=0;
            % ys2=0;
            U0 = U0_(k);
            U1 = U1_(k);
            U0_hbmscs = U0_hbmscs_(k);
            U1_hbmscs = U1_hbmscs_(k);
            U0_far = U0_far_(k);
            theta_far = theta_far_(k);
            D0_x = D0_x_(k);
            D0_y = D0_y_(k);

            % compute the distance in the image for the specific cell
            BW = zeros(size(I0,1),size(I0,2));
            BW(round(max(1,p_cx(k,i-1))),max(1,round(p_cy(k,i-1)))) = 1;
            Dist_n = bwdist(BW); % Gives a matrix with (euclidean) distance from the cell for every pixel
            
            % get position of all cells
            x_position_of_all_cells_with_control = max(1,min(Nx,max(1,round(p_cx(:,i-1)))));% round the x position and assure is between 1 and the image dimension
            y_position_of_all_cells_with_control = max(1,min(Ny,max(1,round(p_cy(:,i-1)))));% round the y position and assure is between 1 and the image dimension
            % first get a matrix with the rows and cols of the distances
            % matrix corresponding to the rows and cols of the cells then
            % take the diagonal which is the actual distance of the single
            % cells. (this way we tackle the need of a for loop)
            dist_nn = diag(Dist_n(x_position_of_all_cells_with_control,y_position_of_all_cells_with_control));
            % in essence we require cells are within a certain interaction distance
            condition_1 = dist_nn<R_int_c;
            % we exclude the cell we are investigating from its own neighbours
            condition_2 = dist_nn > 0;
            % we intersect the two conditions to get the id of actual neighbours of the kth cell
            idx2(k).vec = find(condition_1 & condition_2);
            if isempty(idx2(k).vec)
                xn = 0;
                yn = 0 ;
            end
            
            %% compute index of very far BC cell for each Bc cell

            x_position_of_all_cells_with_control = max(1,min(Nx,max(1,round(p_cx(:,i-1)))));% round the x position and assure is between 1 and the image dimension
            y_position_of_all_cells_with_control = max(1,min(Ny,max(1,round(p_cy(:,i-1)))));% round the y position and assure is between 1 and the image dimension
            % first get a matrix with the rows and cols of the distances
            % matrix corresponding to the rows and cols of the cells then
            % take the diagonal which is the actual distance of the single
            % cells. (this way we tackle the need of a for loop)

            %% HERE insert for loop for the distance of cells outside the image
            % dist_nn = diag(Dist_n(x_position_of_all_cells_with_control,y_position_of_all_cells_with_control));
            % % in essence we require cells are within a certain interaction distance
            % condition_1 = dist_nn>R_int_c_far;
            % % we exclude the cell we are investigating from its own neighbours
            % condition_2 = dist_nn > 0;
            % % we intersect the two conditions to get the id of actual neighbours of the kth cell
            % idx2_far(k).vec = find(condition_1 & condition_2);



            % we want to avoid that there are no neighbours so we increase
            % the interaction radius
            % r=2*R_int_c;
            % while isempty(idx2(k).vec)
            %     condition_1 = dist_nn<r;
            %     idx2(k).vec = find(condition_1 & condition_2);
            %     r=r+1;
            % end

            % if isempty(idx2_far(k).vec)
            %     xn = 0;
            %     yn = 0 ;
            % end
            % compute the squared distance of neighbours cells
            dist_sqr = dist_nn(idx2(k).vec).^2;

            % compute the sum of the reverse of the distances squared
            dens(k)=sum(1./dist_sqr);

            % adjust the volatility according to the local density of the
            % neighbours (between a maximum (sigmax_1, sigmay_1) and a
            % minimum (sigmax_0, sigmay_0) volatility
            sigma_y(k)=min(sigmay_1,max(sigmay_0,sigmay_0+(sigmay_1-sigmay_0)*(1-dens(k)/rho0)));
            sigma_x(k)=min(sigmax_1,max(sigmax_0,sigmax_0+(sigmax_1-sigmax_0)*(1-dens(k)/rho0)));
            
            % HBMSCS
            % compute distance from the kth cell and the hbmscs
            % get position of all hbmscs
            % get position of all cells
            x_position_of_all_hbmscs_cells_with_control = max(1,min(Nx,max(1,round(pos_ix_hbmscs))));% round the x position and assure is between 1 and the image dimension
            y_position_of_all_hbmscs_cells_with_control = max(1,min(Ny,max(1,round(pos_iy_hbmscs))));% round the y position and assure is between 1 and the image dimension
            % first get a matrix with the rows and cols of the distances
            % matrix corresponding to the rows and cols of the cells then
            % take the diagonal which is the actual distance of the single
            % cells. (this way we tackle the need of a for loop)
            dist_nn_hbmscs = diag(Dist_n(x_position_of_all_hbmscs_cells_with_control,y_position_of_all_hbmscs_cells_with_control));
            % in essence we require cells are within a certain interaction distance
            condition_1_hbmscs = dist_nn_hbmscs<R_int_hbmscs;
            % we exclude the cell we are investigating from its own neighbours
            condition_2_hbmscs = dist_nn_hbmscs > 0;
            % we intersect the two conditions to get the id of actual neighbours of the kth cell
            idx2_hbmscs(k).vec = find(condition_1_hbmscs & condition_2_hbmscs);

            % we want to avoid that there are no neighbours so we increase
            % the interaction radius

            % r=2*R_int_hbmscs;
            % while isempty(idx2_hbmscs(k).vec)
            %     condition_1_hbmscs = dist_nn_hbmscs<r;
            %     idx2_hbmscs(k).vec = find(condition_1_hbmscs & condition_2_hbmscs);
            %     r=r+1;
            % end
            if isempty(idx2_hbmscs(k).vec)
                xn_hbmscs = 0;
                yn_hbmscs = 0;
                r_hbmscs_min = a0_hbmscs;
            end
            % compute the squared distance of neighbours cells
            dist_sqr_hbmscs = dist_nn_hbmscs(idx2_hbmscs(k).vec).^2;
            % compute the sum of the reverse of the distances squared
            dens_hbmscs(k)=sum(1./dist_sqr_hbmscs);

            % adjust the volatility according to the local density of the
            % neighbours (between a maximum (sigmax_1, sigmay_1) and a
            % minimum (sigmax_0, sigmay_0) volatility
            sigma_y_hbmscs(k)=min(sigmay_1,max(sigmay_0,sigmay_0+(sigmay_1-sigmay_0)*(1-dens_hbmscs(k)/rho0)));
            sigma_x_hbmscs(k)=min(sigmax_1,max(sigmax_0,sigmax_0+(sigmax_1-sigmax_0)*(1-dens_hbmscs(k)/rho0)));


            % for each neighbour cell of the kth cell compute the xn and yn potential
            for c = 1 : length(idx2(k).vec)
                r_nc(c) = dist_nn(idx2(k).vec(c));
                
                if r_nc(c)< a1 
                    % if distance exceeds a certain threshold both
                    % attractive (linear) and repulsive (exponential) terms
                    % are present
                    GU(c) = -(U0*2 * (r_nc(c)/(a0^2))*exp(-((r_nc(c)/a0)^2)))+U1*2*(a1/r_nc(c));
                else 
                    % if distance is below a certain threshold only
                    % repulsive term is present
                    GU(c) = -U0*2 * (r_nc(c)/(a0^2))*exp(-((r_nc(c)/a0)^2));
                end
                % delta_x(c) = p_cx(idx2(k).vec(c),i-1)-p_cx(k,i-1); % questi sono i delta della cellula i rispetto alle altre vicine
                % delta_y(c) = p_cy(idx2(k).vec(c),i-1)-p_cy(k,i-1);
                % compute relative position of kth cell with respect to its
                % cth neighbour
                delta_x(c) = p_cx(k,i-1)-p_cx(idx2(k).vec(c),i-1); % questi sono i delta della cellula i rispetto alle altre vicine
                delta_y(c) = p_cy(k,i-1)-p_cy(idx2(k).vec(c),i-1);
                alpha2(c) = atan2(delta_y(c),delta_x(c)); % direzione rispetto alla vicina c-esima
                xn(c) = GU(c)*cos(alpha2(c));
                yn(c) = GU(c)*sin(alpha2(c));
                
                % % term to get uniform velocity between neighbouring cells
                % vel_idx2_x=p_cx(idx2(k).vec(c),i-1)-p_cx(idx2(k).vec(c),i-2); % velocità della c-esima in x
                % vel_idx2_y=p_cy(idx2(k).vec(c),i-1)-p_cy(idx2(k).vec(c),i-2); % velocità della c-esima in y
                % xn_unif_v(c)=beta/length(idx2(k).vec)*vel_idx2_x*cos(alpha2(c)); % potenziale di uniformità della vel in x
                % yn_unif_v(c)=beta/length(idx2(k).vec)*vel_idx2_y*sin(alpha2(c)); % potenziale di uniformità della vel in y
            end
            
            clear r_nc delta_x delta_y alpha2 
            % yn_far = U0_far*2 * (p_cy(k,i-1)/(a0_far^2))*exp(-((p_cy(k,i-1)/a0_far)^2));
            
            
            
            % potential far
            % for each neighbour cell of the kth cell compute the xn and yn potential
            % for c = 1 : length(idx2_far(k).vec)
            %     r_nc(c) = dist_nn(idx2_far(k).vec(c));
            % 
            %     % if r_nc(c)< a1_far
            %     % %     if distance exceeds a certain threshold both
            %     % %     attractive (linear) and repulsive (exponential) terms
            %     % %     are present
            %     %     GU(c) = -(U0*2 * (r_nc(c)/(a0^2))*exp(-((r_nc(c)/a0)^2)))+U1*2*(a1/r_nc(c));
            %     % else 
            %         % if distance is below a certain threshold only
            %         % repulsive term is present
            %     GU(c) = +U0_far*2 * (r_nc(c)/(a0_far^2))*exp(-((r_nc(c)/a0_far)^2));
            %     % end
            %     % delta_x(c) = p_cx(idx2(k).vec(c),i-1)-p_cx(k,i-1); % questi sono i delta della cellula i rispetto alle altre vicine
            %     % delta_y(c) = p_cy(idx2(k).vec(c),i-1)-p_cy(k,i-1);
            %     % compute relative position of kth cell with respect to its
            %     % cth neighbour
            %     delta_x(c) = p_cx(k,i-1)-p_cx(idx2_far(k).vec(c),i-1); % questi sono i delta della cellula i rispetto alle altre vicine
            %     delta_y(c) = p_cy(k,i-1)-p_cy(idx2_far(k).vec(c),i-1);
            %     alpha2(c) = atan2(delta_y(c),delta_x(c)); % direzione rispetto alla vicina c-esima
            %     xn_far(c) = GU(c)*cos(alpha2(c));
            %     yn_far(c) = GU(c)*sin(alpha2(c));
            % 
            %     % % term to get uniform velocity between neighbouring cells
            %     % vel_idx2_x=p_cx(idx2(k).vec(c),i-1)-p_cx(idx2(k).vec(c),i-2); % velocità della c-esima in x
            %     % vel_idx2_y=p_cy(idx2(k).vec(c),i-1)-p_cy(idx2(k).vec(c),i-2); % velocità della c-esima in y
            %     % xn_unif_v(c)=beta/length(idx2(k).vec)*vel_idx2_x*cos(alpha2(c)); % potenziale di uniformità della vel in x
            %     % yn_unif_v(c)=beta/length(idx2(k).vec)*vel_idx2_y*sin(alpha2(c)); % potenziale di uniformità della vel in y
            % end
            % for each neighbour cell of the kth cell compute the xn and yn potential
            for c = 1 : length(idx2_hbmscs(k).vec)
                r_nc(c) = dist_nn_hbmscs(idx2_hbmscs(k).vec(c));
                
                if r_nc(c)< a0_hbmscs 
                    % if distance exceeds a certain threshold both
                    % attractive (linear) and repulsive (exponential) terms
                    % are present
                    U0_far = U0_far_(k)/10;
                    GU_hbmscs(c) = -(U1_hbmscs*2 * (r_nc(c)/(a0_hbmscs^2))*exp(-((r_nc(c)/a0_hbmscs)^2)))+U0_hbmscs*2*((r_nc(c)/a1_hbmscs)^(-1));
                else 
                    % if distance is below a certain threshold only
                    % repulsive term is present
                    
                    GU_hbmscs(c) = -U1_hbmscs*2 * (r_nc(c)/(a0_hbmscs^2))*exp(-((r_nc(c)/a0_hbmscs)^2));
                end
                
                % delta_x(c) = p_cx(idx2(k).vec(c),i-1)-p_cx(k,i-1); % questi sono i delta della cellula i rispetto alle altre vicine
                % delta_y(c) = p_cy(idx2(k).vec(c),i-1)-p_cy(k,i-1);
                % compute relative position of kth cell with respect to its
                % cth neighbour
                delta_x_hbmscs(c) = p_cx(k,i-1)-pos_ix_hbmscs(idx2_hbmscs(k).vec(c)); % questi sono i delta della cellula i rispetto alle altre vicine
                delta_y_hbmscs(c) = p_cy(k,i-1)-pos_ix_hbmscs(idx2_hbmscs(k).vec(c));
                alpha2_hbmscs(c) = atan2(delta_y_hbmscs(c),delta_x_hbmscs(c)); % direzione rispetto alla vicina c-esima
                xn_hbmscs(c) = GU_hbmscs(c)*cos(alpha2_hbmscs(c));
                yn_hbmscs(c) = GU_hbmscs(c)*sin(alpha2_hbmscs(c));
                r_hbmscs_min = min(r_nc);
                % % term to get uniform velocity between neighbouring cells
                % vel_idx2_x=p_cx(idx2(k).vec(c),i-1)-p_cx(idx2(k).vec(c),i-2); % velocità della c-esima in x
                % vel_idx2_y=p_cy(idx2(k).vec(c),i-1)-p_cy(idx2(k).vec(c),i-2); % velocità della c-esima in y
                % xn_unif_v(c)=beta/length(idx2(k).vec)*vel_idx2_x*cos(alpha2(c)); % potenziale di uniformità della vel in x
                % yn_unif_v(c)=beta/length(idx2(k).vec)*vel_idx2_y*sin(alpha2(c)); % potenziale di uniformità della vel in y
            end
            
            factor = min(r_hbmscs_min/a0_hbmscs,1);
            clear r_nc alpha2_hbmscs delta_x_hbmscs delta_y_hbmscs
            
            % yn_far = U0_far*2 * (abs((2*Ny-p_cy(k,i-1)))/2*Ny)*factor;
            xn_far = U0_far*2*(1-exp(-((2*Ny-p_cy(k,i-1)/10)^2)))*factor*sin(theta_far);
            yn_far = U0_far*2*(1-exp(-((2*Ny-p_cy(k,i-1)/10)^2)))*factor*cos(theta_far);
            % sum all the potentials (their gradients
            pot_xn = sum(xn) + sum(xn_hbmscs) + xn_far;%+ sum(xn_far);
            pot_yn = sum(yn) + sum(yn_hbmscs) + yn_far;
                
            vx_random = D0_x*(1+0.2*randn());
            vy_random = D0_y*(1+0.2*randn());

            % pot_xn = sum(xn_hbmscs);
            % pot_yn = sum(yn_hbmscs);
            % pot_xn = 0;
            % pot_yn = 0;

            % apply adams method to update vx and vy
            vsol_x_cand = vel_x(k,i-1) + pot_xn + D0_x*(randn());
            vsol_y_cand = vel_y(k,i-1) + pot_yn + D0_y*(randn());

            vsol_x = (1-alpha)*vsol_x_cand;
            vsol_y = (1-alpha)*vsol_y_cand;
            % vsol_x = vel_x(k,i-1) - alpha*vel_x(k,i-1) + pot_xn;
            % vsol_y = vel_y(k,i-1) - alpha*vel_y(k,i-1) + pot_yn;

            % [ t_x , vsol_x ] = Adams_esplicito ( @( t , v) (-alpha*v+pot_xn) , vel_x(k,i-1), [i-1 i], 1 );
            % 
            % [ t_y , vsol_y ] = Adams_esplicito ( @( t , v) (-alpha*v+pot_yn) , vel_y(k,i-1), [i-1 i], 1 );
            
            % insert the velocity as delta in the position
            % vel_x(k,i) = vsol_x(end);
            % vel_y(k,i) = vsol_y(end);
            dx(k) = vsol_x(end);
            dy(k) = vsol_y(end);

            clear delta_x delta_y delta_x_hbmscs delta_y_hbmscs alpha2 alpha2_hbmscs GU GU_hbmscs BW Dist_n bw r_nc GS alphas xs ys;
            % p_sp(k).y(i) = min(Ny-2*cell_dim,dy(k)+p_sp(k).y(i-1)+D0+(sigma_y(k)+sigma_y_hbmscs(k))*randn());
            % p_sp(k).x(i) = min(Nx-2*cell_dim,dx(k)+p_sp(k).x(i-1)+D0+(sigma_x(k)+sigma_x_hbmscs(k))*randn());
            
            % p_sp(k).y(i) = min(Ny-2*cell_dim,dy(k)+p_sp(k).y(i-1)+D0+(sigma_y(k))*randn());
            % p_sp(k).x(i) = min(Nx-2*cell_dim,dx(k)+p_sp(k).x(i-1)+D0+(sigma_x(k))*randn());
            p_sp(k).y(i) = min(Ny-2*cell_dim,dy(k)+p_sp(k).y(i-1)+D0_y*0.2*(randn()));%+D0_y*(randn()));%+2*rand()));
            p_sp(k).x(i) = min(Nx-2*cell_dim,dx(k)+p_sp(k).x(i-1)+D0_x*0.2*(randn()));%+D0_x*(randn()));
            p_sp(k).U0_far_x(i) = U0_far*sin(theta_far); 
            p_sp(k).U0_far_y(i) = U0_far*cos(theta_far);
            p_sp(k).U0_far(i) = U0_far;
            p_sp(k).theta_far(i) = theta_far;
            p_sp(k).vx_random(i) = vx_random;
            p_sp(k).vy_random(i) = vy_random;
            p_sp(k).U0(i) = U0;
            p_sp(k).U1(i) = U1;
            p_sp(k).U0_hbmscs(i) = U0_hbmscs;
            p_sp(k).U1_hbmscs(i) = U1_hbmscs;


            vel_x(k,i) = p_sp(k).x(i)-p_sp(k).x(i-1);
            vel_y(k,i) = p_sp(k).y(i)-p_sp(k).y(i-1);
            
        
            p_sp(k).t(i) = i;

            

        end
        if (visualize)&&(mod(i,1)==0)

                imshow(I0(:,:));

                for n_cell = 1 : length(p_sp)
                    hold on
            
                    viscircles([p_sp(n_cell).y(i),p_sp(n_cell).x(i)],20,'Color','g');
                    % plot(p_sp(n_cell).y(t),p_sp(n_cell).x(t),'o','Color','g','MarkerSize',10)
                end
                for n_cell = 1 : length(pos_ix_hbmscs)
                    hold on
                    viscircles([pos_iy_hbmscs(n_cell),pos_ix_hbmscs(n_cell)],30);
                    % viscircles(pos_iy_hbmscs(n_cell),pos_ix_hbmscs(n_cell)],20)
                    % plot(pos_iy_hbmscs(n_cell),pos_ix_hbmscs(n_cell),'o','MarkerSize',30,'Color','r');
                end
                title("Time: "+num2str(i))
                drawnow
        end
    end
    clear k 
    save(name_save)
    % save(name_save,'p_sp','I0','pos_ix_hbmscs','pos_iy_hbmscs','U0_','U1_',"U0_hbmscs_","U1_hbmscs_","-v7.3")
end





        




