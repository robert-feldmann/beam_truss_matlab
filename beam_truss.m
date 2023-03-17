classdef beam_truss < handle

    properties
        P           % truss node matrix (main nodes)        
        Beam          % truss beam meatrix
        K           % stiffness matrix components
        M           % mass matrix components
        q_glob      % vector of generalized coordinates / node information matrix 
        param       % truss parameters 
        n           % number of elements per beam
        koinz       % coincidence matrix
        node_coord  % coordinate matrix of FE nodes

        V_norm      % matrix of normalized modes
        f_eig       % vector of eigenfrequencies 

    end

    methods
        
        function obj = beam_truss(P,Beam,n)
            % constructor method for instance of class beam truss

            % inputs: 
            % matrix of truss nodes (see readme)
            % matrix of truss beams (see readme)
            % number of elements per beam
            
            % inputs
            obj.P = P;
            obj.Beam = Beam;
            obj.n = n;

            % parameters (to be adapted accordingly
            obj.param.beam.d =  0.005*2;                                        % beam diameter
            obj.param.beam.E = 71.8e9;                                          % Young's modulus in N/m^2 (Pa) 
            obj.param.beam.rho = 2813;                                          % density aluminum in kg/m^3 
            obj.param.ball.r = 0.03;                                            % connecting ball radius 
            obj.param.ball.m = 0.268;                                           % connecting ball mass 
            obj.param.beam.G = obj.param.beam.E/(2*(1+0.34));                   % shear modulus in N/m^2
            obj.param.beam.A = (pi*obj.param.beam.d^2)/4;                       % cross sectional area in m^2
            obj.param.beam.Iy = pi/4*(obj.param.beam.d/2)^4;                    % moment of inertia in m^4 
            obj.param.beam.Iz = obj.param.beam.Iy;                              % moment of inertia in m^4
            obj.param.beam.It = obj.param.beam.Iy + obj.param.beam.Iz;          % polar moment of inertia in m^4
            obj.param.beam.theta = obj.param.beam.A * obj.param.beam.rho*obj.param.beam.d^2/8; % mass moment of inertia/length in kg*m
            obj.param.beam.m = obj.param.beam.rho*obj.param.beam.A;             % mass matrix parameters
            obj.param.ball.theta = 0.4* obj.param.ball.m* obj.param.ball.r^2;   % mass moment of inertia/length in kg*m (assuming full ball)

           
            % creation of the system matrices 
            setup_beam_truss(obj);
            
        end

        function setup_beam_truss(obj)
            
            % The derivation of the mass and stiffness matrix in this method is based on the books

                % Géradin, Michel, and Daniel J. Rixen. Mechanical vibrations: theory and application to structural dynamics. John Wiley & Sons, 2014.
                % Shabana, Ahmed A. Dynamics of multibody systems. Cambridge university press, 2020.

            dof = 6; % number of dof: 3 translational, 3 rotational
            
            %% -1- generating the element matrices
            disp('generating element matrices...')
            % shape matrice 
            syms xsy lsy 
            ksi=xsy/lsy;
            F=[1-ksi; ksi; 1-3*ksi^2+2*ksi^3 ;lsy*ksi*(1-ksi)^2; ksi^2*(3-2*ksi); lsy*ksi^2*(ksi-1)];
            S=[F(1)                               0                             0                  0;
                0                                  F(3)                          0                  0;
                0                                  0                             F(3)               0;
                0                                  0                             0                  F(1);
                0                                  0                             -F(4)              0;
                0                                  F(4)                          0                  0;
                F(2)                               0                             0                  0;
                0                                  F(5)                          0                  0;
                0                                  0                             F(5)               0;
                0                                  0                             0                  F(2);
                0                                  0                             -F(6)              0;
                0                                  F(6)                          0                  0].';

            % calculate element mass and stiffness matrices

            % element mass matrix: translation
            tmp1=blkdiag(1,1,1,0);
            tmp2=S.'*tmp1*S;
            Me_t=eval(int(tmp2,xsy,0,lsy));
            clear tmp1 tmp2

            % element mass matrix: rotation
            tmp1=blkdiag(0,0,0,1);
            tmp2=S.'*tmp1*S;
            Me_r=eval(int(tmp2,xsy,0,lsy));
            clear tmp1 tmp2

            DS= [diff(S(1,:),xsy); % differentiation of form function
                diff(S(2,:),xsy,2);
                diff(S(3,:),xsy,2);
                diff(S(4,:),xsy,1)];

            % element stiffness matrix: stress
            tmp1=blkdiag(1, 0, 0, 0);
            tmp2=DS.'*tmp1*DS;
            Ke_zd=eval(int(tmp2,xsy,0,lsy));
            clear tmp1 tmp2

            % element stiffness matrix: bending
            tmp1=blkdiag(0, 1, 1, 0);
            tmp2=DS.'*tmp1*DS;
            Ke_b=eval(int(tmp2,xsy,0,lsy));
            clear tmp1 tmp2

            % element stiffness matrix: torsion
            tmp1=blkdiag(0, 0, 0, 1);
            tmp2=DS.'*tmp1*DS;
            Ke_to=eval(int(tmp2,xsy,0,lsy));
            clear tmp1 tmp2

            % additional translational and rotational inertia of connecting balls
            M_ball_t = blkdiag(ones(3,3),zeros(3,3));
            M_ball_r = blkdiag(zeros(3,3),ones(3,3));

            clear xsy lsy Esy Gsy Asy Isy Itsy msy thetasy ksi rsy
            clear Ke S DS F

            % note: The lengths are not inserted in this section because the beams have different lengths in some cases. They are only in section 5.

            %% -2- diskretization
            disp('discretization...')
            % K struct:
                % bar number (see matrix Beam)
                % element node type:
                % 1 - main node
                % ... more types of node types can be implemented and inserted here, for example connecting elements between the beams and the connecting balls 
                % 4 - beam B
                % number (consecutive)
                % coordinates

            % node_index: struct that contains the node indices of the bars.
            % bar number (see matrix Beam),
            % element type,
            % element number
            
            [obj.node_coord,node_index] =  beam_truss.truss_discretization(obj.P,obj.Beam,obj.n);
            
            %% -3- coincidence matrices (referring to beam nodes)
            disp('generating coincidence matrices and index matrix...')
            % In the following, coincidence matrices are created specifically for each beam
            % Each line describes the nodes of a finite element with respect to beam number, node type and node number 

            koinz_var_names = {'beam_start', 'node_type_start', 'node_number_start','beam_end', 'beam_type_end', 'node_number_end'};

            for bl = 1:size(obj.Beam,1)

                % coincidence matrix for beams
                obj.koinz.B.Balken{bl} = ...
                    array2table(...
                    [table2array(node_index.beam{bl}(1:end-1,:)), table2array(node_index.beam{bl}(2:end,:))], ... % each node is coupled with the respective following node.
                    'VariableNames', koinz_var_names);

            end

            %% -4- index matrix for vector of generalized coordinates for each beam (q_loc_index)
            
            % structure
            % beam number (see matrix Beam),
            % element type (1 - main nodes, 4 - beam B )
            % element number
            % dof number: 1 - x,  2 - y, 3 - z, 4 - \psi_x, 5 - \psi_y, 6 - \psi_z

            for bl = 1: size(obj.Beam,1) % loop over beams

                n_kno = size(node_index.beam{bl},1); % number of nodes for beam 
                q_loc_index.Balken{bl} = ...
                    [reshape(repmat(table2array(node_index.beam{bl}(:,1)).',dof,1),[],1),...
                    reshape(repmat(table2array(node_index.beam{bl}(:,2)).',dof,1),[],1),...
                    reshape(repmat(table2array(node_index.beam{bl}(:,3)).',dof,1),[],1),...
                    repmat(linspace(1,dof,dof).',n_kno,1 )];

                % table 
                q_loc_index.Balken{bl} = array2table( q_loc_index.Balken{bl}, ...
                    'VariableNames', {'beam_number', 'node_type', 'node_number', 'dof_number'});

                % meta data
                q_loc_index.Balken{bl}.Properties.VariableDescriptions{'dof_number'} = ...
                    '1 - x,  2 - y, 3 - z, 4 - \psi_x, 5 - \psi_y, 6 - \psi_z';
                q_loc_index.Balken{bl}.Properties.VariableDescriptions{'node_number'} = ...
                    'consecutive number of the element of the respective element type';
                q_loc_index.Balken{bl}.Properties.VariableDescriptions{'node_type'} = ....
                    '1 - main node, 2 - beam node B';

            end
            clear bl

            %% -5- index matrix for global vector of generalized coordinates

            % in the following, the index matrix q_glob is created. It is the extension of q_loc_index (specific for each beam).
            % q_glob refers to the FG of the whole structure.

            % structure index matrix q_glob
            % bar number (see matrix Beam and ) %Note: main nodes have bar number 0
            % element type, (as above)
            % element number,
            % dof number: 1 - x, 2 - y, 3 - z, 4 - psi_x, 5 - psi_y, 6 - psi_z

            % structure: 1. dof of the main nodes
            %            2. dof of the beams in the order of the matrix Beam

            % assumption: each beam has two main nodes (matrix Beam).

            n_kno_HK = size(obj.P,1); % number of main nodes 
            % 1. dof of the main nodes
            obj.q_glob = [   reshape(repmat(zeros(n_kno_HK,1),dof,1),[],1),... (beam number (main node): 0)
                reshape(repmat(ones(n_kno_HK,1).',dof,1),[],1),... (element type: main node, 1)
                reshape(repmat(obj.P(:,1).',dof,1),[],1),... (element numbers)
                repmat(linspace(1,dof,dof).', n_kno_HK, 1 )  ]; 
            
            % 2. dof of the beams in the order of the matrix Beam
            % q_glob is successively extended by bar dofs in the loop
            for bl = 1: size(obj.Beam,1) % Schleife über alle Balken
                obj.q_glob = [ obj.q_glob ; ...
                    table2array(q_loc_index.Balken{bl}(dof+1:end-dof,:))]; % start and end nodes (=main nodes) are cleared out
            end

            obj.q_glob = array2table(obj.q_glob, 'Variablenames', {'beam_number', 'node_type', 'node_number', 'dof_number'});
            clear n_kno_HK bl

            %% -6- generation of the mass ans stiffness matrices
            disp('assembly of mass and stiffness matrices...')

            % The matrices M_ff and K_ff describe the inertia and stiffness of the structure in the
            % solid coordinate system with respect to the nodal coordinates.
            
            % initialization
            M_ff=zeros(size(obj.q_glob,1));
            K_ff=zeros(size(obj.q_glob,1));
 
            % stiffness matrix
            K_ff_Bzd=K_ff; 
            K_ff_Bb=K_ff;
            K_ff_Bto=K_ff;
            
            % mass matrix
            M_ff_Bt=M_ff;
            M_ff_Br=M_ff;
            M_ff_ball_t=M_ff;
            M_ff_ball_r=M_ff;

            clear M_ff K_ff M_rf
            w = waitbar(0,'assembly in progress...'); 
            
            for bl=1:size(obj.Beam,1)
                
                %transformation of the element deformation (u,v,w) into the structural coordinate system 
                C1.Balken{bl}=zeros(3);
                C2.Balken{bl}=zeros(12);
                C2.Balken{bl}=inv(beam_truss.Trans_loc_glob(obj.P(obj.Beam(bl,2),2:4),obj.P(obj.Beam(bl,3),2:4),dof));
                C1.Balken{bl}=inv(C2.Balken{bl}(1:3,1:3)); % todo: Trafo vereinfachen 

                for el_B = 1:size(obj.koinz.B.Balken{bl},1) % loop over beam elements 

                    % matrix B1 maps the element dof to their respective position in q_glob
                    B1=zeros(12,size(obj.q_glob,1));
                    
                    % if query because of indexing convention
                    if table2array(obj.koinz.B.Balken{bl}(el_B,2)) == 1 
                        [~,tmp1,~]=beam_truss.intersectrows_fast(table2array(obj.q_glob),([0 1 obj.P(obj.Beam(bl,2)) 1]));
                    else 
                        [~,tmp1,~]=beam_truss.intersectrows_fast(table2array(obj.q_glob),([table2array(obj.koinz.B.Balken{bl}(el_B,1:3)),1]));
                    end
                    if table2array(obj.koinz.B.Balken{bl}(el_B,5)) == 1
                        [~,tmp2,~]=beam_truss.intersectrows_fast(table2array(obj.q_glob),([0 1 obj.P(obj.Beam(bl,3)) 1]));
                    else
                        [~,tmp2,~]=beam_truss.intersectrows_fast(table2array(obj.q_glob),([table2array(obj.koinz.B.Balken{bl}(el_B,4:6)),1]));
                    end
                    B1(1:6,tmp1:tmp1+5)=blkdiag(1,1,1,1,1,1);
                    B1(7:12,tmp2:tmp2+5)=blkdiag(1,1,1,1,1,1);

                    % assembly of mass ans stiffness matrices using loops
                    % over beams and elements

                    % mass matrix: torsion part
                    M_B_t = beam_truss.M_Matrix(obj.node_coord.Balken{bl}.el_type{ table2array( obj.koinz.B.Balken{bl}(el_B,2)) }.num{ table2array(obj.koinz.B.Balken{bl}(el_B,3)) }.coordinates,... % Anfangsknoten
                        obj.node_coord.Balken{bl}.el_type{ table2array( obj.koinz.B.Balken{bl}(el_B,5)) }.num{ table2array(obj.koinz.B.Balken{bl}(el_B,6)) }.coordinates,... % Endknoten
                        Me_t,C2.Balken{bl},B1);
                    M_ff_Bt=M_ff_Bt+M_B_t;
                    
                    % mass matrix rotation part
                    M_B_r = beam_truss.M_Matrix(obj.node_coord.Balken{bl}.el_type{ table2array( obj.koinz.B.Balken{bl}(el_B,2)) }.num{ table2array(obj.koinz.B.Balken{bl}(el_B,3)) }.coordinates,... % Anfangsknoten
                        obj.node_coord.Balken{bl}.el_type{ table2array( obj.koinz.B.Balken{bl}(el_B,5)) }.num{ table2array(obj.koinz.B.Balken{bl}(el_B,6)) }.coordinates,... % Endknoten
                        Me_r,C2.Balken{bl},B1);
                    M_ff_Br=M_ff_Br+M_B_r;

                    % stiffness matrix stress part
                    K_B_zd = beam_truss.M_Matrix(obj.node_coord.Balken{bl}.el_type{ table2array( obj.koinz.B.Balken{bl}(el_B,2)) }.num{ table2array(obj.koinz.B.Balken{bl}(el_B,3)) }.coordinates,... % Anfangsknoten
                        obj.node_coord.Balken{bl}.el_type{ table2array( obj.koinz.B.Balken{bl}(el_B,5)) }.num{ table2array(obj.koinz.B.Balken{bl}(el_B,6)) }.coordinates,... % Endknoten
                        Ke_zd,C2.Balken{bl},B1);
                    K_ff_Bzd=K_ff_Bzd+K_B_zd;
                    
                    % stiffness matrix bending part
                    K_B_b = beam_truss.M_Matrix(obj.node_coord.Balken{bl}.el_type{ table2array( obj.koinz.B.Balken{bl}(el_B,2)) }.num{ table2array(obj.koinz.B.Balken{bl}(el_B,3)) }.coordinates,... % Anfangsknoten
                        obj.node_coord.Balken{bl}.el_type{ table2array( obj.koinz.B.Balken{bl}(el_B,5)) }.num{ table2array(obj.koinz.B.Balken{bl}(el_B,6)) }.coordinates,... % Endknoten
                        Ke_b,C2.Balken{bl},B1);
                    K_ff_Bb=K_ff_Bb+K_B_b;
                    
                    % stiffness matrix torsion part
                    K_B_to = beam_truss.M_Matrix(obj.node_coord.Balken{bl}.el_type{ table2array( obj.koinz.B.Balken{bl}(el_B,2)) }.num{ table2array(obj.koinz.B.Balken{bl}(el_B,3)) }.coordinates,... % Anfangsknoten
                        obj.node_coord.Balken{bl}.el_type{ table2array( obj.koinz.B.Balken{bl}(el_B,5)) }.num{ table2array(obj.koinz.B.Balken{bl}(el_B,6)) }.coordinates,... % Endknoten
                        Ke_to,C2.Balken{bl},B1);
                    K_ff_Bto=K_ff_Bto+K_B_to;

                end
                waitbar( bl / size(obj.Beam,1), w)
            end
            
            close(w)
            
            % mass supplement for sphere
            % matrix M_Sphere has only diagonal elements since all axes are always principal axes
    
            for HK = 1 : size(obj.P,1)  % loop over main nodes
                for i_ball = 1:6 % i_kugel = diagonal element of the matrix M_sphere
                    [~, i_glob ,~] =  intersect( ...
                        table2array(obj.q_glob) , ...
                        [0, 1, HK, i_ball],... %beam number 0 (main node), node type 1 (main node), node number , dof number)
                        'rows');
                    M_ff_ball_t(i_glob, i_glob) = M_ff_ball_t(i_glob, i_glob) + M_ball_t(i_ball, i_ball);
                    M_ff_ball_r(i_glob, i_glob) = M_ff_ball_r(i_glob, i_glob) + M_ball_r(i_ball, i_ball);

                end
            end

            %% merge partial matrices into struct

            obj.K.Bzd = sparse(K_ff_Bzd); % sparse matrices speed up calculation
            obj.K.Bb = sparse(K_ff_Bb);
            obj.K.Bto = sparse(K_ff_Bto);
            obj.M.ball_t = sparse(M_ff_ball_t);
            obj.M.Bt = sparse(M_ff_Bt);
            obj.M.ball_r = sparse(M_ff_ball_r);
            obj.M.Br = sparse(M_ff_Br);

            obj.K.K_total = obj.param.beam.E*obj.param.beam.A*obj.K.Bzd + obj.param.beam.E*obj.param.beam.Iz*obj.K.Bb + obj.param.beam.G*obj.param.beam.It*obj.K.Bto;
            
            obj.M.M_total = obj.param.beam.m *obj.M.Bt + obj.param.ball.m*obj.M.ball_t + obj.param.beam.theta *obj.M.Br + obj.param.ball.theta *obj.M.ball_r;
            
            disp('***system generation complete***')
            
            plot_geom(obj)

        end

        function plot_geom(obj)

            disp('plotting geometry...')

            % function to plot the geometry     

            close all; 
            figure; hold on;
            
            w = waitbar(0,'plotting geometry...');

            % plot beam elements
            for i_b = 1:size(obj.Beam, 1) % loop over beams
                for i_el = 1:size(obj.koinz.B.Balken{1,i_b},1) % loop over elements
                    id_anf = table2array(obj.koinz.B.Balken{1,i_b}(i_el,1:3));
                    id_end = table2array(obj.koinz.B.Balken{1,i_b}(i_el,4:6));
                    p = [ obj.node_coord.Balken{1,id_anf(1)}.el_type{1,id_anf(2)}.num{1,id_anf(3)}.coordinates;  ...
                        obj.node_coord.Balken{1,id_end(1)}.el_type{1,id_end(2)}.num{1,id_end(3)}.coordinates];
                    scatter3( p(:,1), p(:,2), p(:,3) ,'k', 'filled', 'LineWidth', 0.01); % plot beam nodes (todo: 
                    line( p(:,1), p(:,2), p(:,3) , 'LineWidth', 2, 'Color','green') % plot beams
                end
                drawnow
                waitbar( i_b / size(obj.Beam,1), w)
            end
            
            close(w)

            % main nodes: bold
            scatter3( obj.P(:,2), obj.P(:,3), obj.P(:,4),'r', 'filled', 'LineWidth', 2);
            shg

            ax = gca;
            ax.Visible = false;
            view([0 1 0])  
            axis equal

            % numbering of main nodes 
            text( obj.P(:,2), obj.P(:,3), obj.P(:,4), num2str(obj.P(:,1)), 'FontSize',10, 'Color','red');
            
            % numbering of beams
            for bl = 1: size(obj.Beam,1)
                text_loc = [beam_truss.Kn_ele( [obj.P(obj.Beam(bl,2),2:end), obj.P(obj.Beam(bl,3), 2:end), 2])];
                % note: Kn_ele is a function that is used here to get the text to the center of the bar.
                text( text_loc(1), text_loc(2), text_loc(3), strcat(num2str(bl)), 'FontSize',10 , 'Color', 'blue' );
            end

        end

        function modal_analysis(obj)
            % conduct a modal analysis

            disp('conducting modal analysis...')

            % make matrices symmetrical (numerical error induced in trans_loc_glob)
            K_glob_T = (obj.K.K_total + obj.K.K_total')/2;
            M_glob_T = (obj.M.M_total + obj.M.M_total')/2;
            
            num_eig = 35; % the number of eigenvectors to be calculated
            [V,D] = eigs(K_glob_T, M_glob_T, num_eig , 0); 
            
            % normalization
            obj.V_norm=zeros(size(V,1), size(D,1)-6); %rigid body modes are not examined
            for idx = 7: size(D,1)
                obj.V_norm(:,idx) = ( V(:,idx).'*M_glob_T*V(:,idx) ) \ V(:,idx);
            end
            clear idx;
                
            % eigenfrequencies 
            obj.f_eig = sqrt(diag(D))/(2*pi); 

        end
        
        function disp_modal_analysis(obj, i_mode)
            % visualize the calculated modes

            disp('visualization eigenmodes');

             % selection eigenmode (7 - first elastic eigenmode)
            i_el_mode = i_mode - 6;
            knot_size = 0.005; % display size of the nodes
            n_frames = 30; % number of frames 

            % init video
            scaleFact = 5e-1/max(obj.V_norm(:,i_mode)); % scaling of the eigenmodes
            t=linspace(0, 1, n_frames); scale=scaleFact*sin(t*2*pi);

            video_name = strcat(string(datetime('now','Format','yymmdd_HHmm_')), 'mode_', num2str(i_el_mode) );
            v = VideoWriter(video_name,...
                'Motion JPEG AVI'); % Note: alter video format here
            v.FrameRate = n_frames/3;
            %     v.LosslessCompression = 'true';
            %     v.Quality =100;
            open(v);

            w = waitbar(0,'Creating the frames for modes visualization...');

            for idx=1:n_frames

                %% plot settings
                waitbar( idx / n_frames, w)

                fig = figure(idx); hold on
                fig.Units = 'normalized';
                fig.OuterPosition = [0 0 1 1];

                % scaling of the eigenmodes
                V = scale(idx)*obj.V_norm(:,i_mode);

                %% plot

                % beam elements
                for i_b = 1:size(obj.Beam, 1) % loop over beams
                    for i_el = 1:size(obj.koinz.B.Balken{1,i_b},1) % loop over elements
                        
                        if table2array(obj.koinz.B.Balken{1,i_b}(i_el,2)) == 1 % if query due to indexing convention 
                            id_anf = [0 1 obj.P(obj.Beam(i_b,2))];
                            p_anf = obj.P(obj.Beam(i_b,2),2:end);
                            scatter3(p_anf(1), p_anf(2), p_anf(3), 'k', 'filled', 'LineWidth', 0.01); % plot main nodes
                        else
                            id_anf = table2array(obj.koinz.B.Balken{1,i_b}(i_el,1:3));
                            p_anf = obj.node_coord.Balken{1,id_anf(1)}.el_type{1,id_anf(2)}.num{1,id_anf(3)}.coordinates; 
                        end
                        if table2array(obj.koinz.B.Balken{1,i_b}(i_el,5)) == 1
                            id_end = [0 1 obj.P(obj.Beam(i_b,3))];
                            p_end = obj.P(obj.Beam(i_b,3),2:end);  
                            scatter3(p_end(1), p_end(2), p_end(3), 'k', 'filled', 'LineWidth', 0.01); % plot main nodes
                        else
                            id_end = table2array(obj.koinz.B.Balken{1,i_b}(i_el,4:6));
                            p_end = obj.node_coord.Balken{1,id_end(1)}.el_type{1,id_end(2)}.num{1,id_end(3)}.coordinates; 
                        end
                        
                        p = [p_anf; p_end];
                        
                        % Plot displacements
                        % x - coordinate (element start)
                        [~, i_row,~] =  beam_truss.intersectrows_fast( table2array(obj.q_glob), ...
                            [ id_anf(1), id_anf(2), id_anf(3), 1], 'rows');
                        p(1,1) = p(1,1) + V(i_row);
                        % y - coordinate (element start)
                        [~, i_row,~] =  beam_truss.intersectrows_fast( table2array(obj.q_glob), ...
                            [ id_anf(1), id_anf(2), id_anf(3), 2], 'rows');
                        p(1,2) = p(1,2) + V(i_row);
                        % z - coordinate (element start)
                        [~, i_row,~] =  beam_truss.intersectrows_fast( table2array(obj.q_glob), ...
                            [ id_anf(1), id_anf(2), id_anf(3), 3], 'rows');
                        p(1,3) = p(1,3) + V(i_row);
                        % x - coordinate (element end)
                        [~, i_row,~] =  beam_truss.intersectrows_fast( table2array(obj.q_glob), ...
                            [ id_end(1), id_end(2), id_end(3), 1], 'rows');
                        p(2,1) = p(2,1) + V(i_row);
                        % y - coordinate (element end)
                        [~, i_row,~] =  beam_truss.intersectrows_fast( table2array(obj.q_glob), ...
                            [ id_end(1), id_end(2), id_end(3), 2], 'rows');
                        p(2,2) = p(2,2) + V(i_row);
                        % z - coordinate (element end)
                        [~, i_row,~] =  beam_truss.intersectrows_fast( table2array(obj.q_glob), ...
                            [ id_end(1), id_end(2), id_end(3), 3], 'rows');
                        p(2,3) = p(2,3) + V(i_row);

                        scatter3( p(:,1), p(:,2), p(:,3) ,'k', 'filled', 'LineWidth', knot_size); % todo
                        line( p(:,1), p(:,2), p(:,3) , 'LineWidth', 2, 'Color','green')  
                    end
                end
                
                %% plot settings

                % xlim([-0.75 0.75]); % uncomment and adapt if necessary
                % ylim([-0.75 0.75]);
                % zlim([ 0.3 1]);
                % axis equal 
                axis off
                view([0 0 1]); % important: here you can alter the view of the truss
                text(-0.75, -0.75, 0.3 , strcat('eigenfrequency $$f_{',num2str(i_mode),"}$$= " , num2str(obj.f_eig(i_mode)), " Hz" ), 'Fontsize', 20, 'Interpreter','latex');

                %% write video

                frame = getframe(gcf);
                writeVideo(v,frame);
                close(fig); clear fig; 

            end

            close(w); % waitbar
            close(v); % save video

            implay(strcat(video_name, '.avi'));

        end
        
    end

    methods(Static)

        function T_ges = Trans_loc_glob( kn_a, kn_e, dof ) 
            % The function derives the tansformation matrix from the local 
            % coordinate system of the beam into the solid coordinate
            % system. 

            % kn_a, kn_e are the position vectors of the nodes of an element

            % The local system(e_x,e_y,e_z) is defined so that e_x is the longitudinal direction of the beam and (e_x,e_y,e_z)
            % form a right system

            % The structural coordinate system is defined by (e_X,e_Y,e_Z) with e_X=(1 0 0).' e_Y=(0 1 0).' e_Z=(0 0 1).'

            x1 = kn_a(1);
            y1 = kn_a(2);
            z1 = kn_a(3);

            x2 = kn_e(1);
            y2 = kn_e(2);
            z2 = kn_e(3);

            % element length
            le = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);

            %% element basis vectors in structural coordinate system

            % assumption: bar is round, only orientation x-axis important

            e_x = 1/le.*[x2-x1; y2-y1; z2-z1];

            % (local y-axis should point into global z-axis)

            if cond([e_x,[0 0 1]']) < cond([e_x,[0 1 0]'])
                e_y_tmp = [0 0 1]';
            else
                e_y_tmp = [0 1 0]'; % applies in certain cases of beam orientations
            end

            e_z = cross(e_x, e_y_tmp); e_z = e_z./norm(e_z,2);
            e_y = -cross(e_x, e_z); e_y = e_y./norm(e_y,2);

            %% basis vectors structural coordinate system

            e_X = [ 1 0 0]';
            e_Y = [ 0 1 0]';
            e_Z = [ 0 0 1]';

            %% transformation matrix (see Rixen)
            
            T(1,1) = dot(e_x, e_X); 
            T(2,1) = dot(e_x, e_Y);
            T(3,1) = dot(e_x, e_Z);
            
            T(1,2) = dot(e_y, e_X);
            T(2,2) = dot(e_y, e_Y);
            T(3,2) = dot(e_y, e_Z);
            
            T(1,3) = dot(e_z, e_X);
            T(2,3) = dot(e_z, e_Y);
            T(3,3) = dot(e_z, e_Z);
                        
            if dof == 6
                T_ges = [ T        zeros(3) zeros(3) zeros(3);...
                    zeros(3) T        zeros(3) zeros(3);...
                    zeros(3) zeros(3) T        zeros(3);...
                    zeros(3) zeros(3) zeros(3) T        ];
            elseif dof == 3
                T_ges = [ T        zeros(3) ;...
                    zeros(3) T        ];
            end
        end

        function [K,node_index] = truss_discretization(P,S,n)
            % function for discretization of the beams
            % output is a struct with the coordinates and index matrix 

            node_index_table_names = {'beam', 'element_type', 'node_number'};

            for b1 = 1:size(S,1) % loop over beams
                % main nodes (start)
                K.Balken{b1}.el_type{1}.num{1}.coordinates = P(S(b1,2),2:4);
                node_index.beam{b1} = [ ];
                node_index.beam{b1} = array2table([b1 , 1 , 1 ], ...
                    'VariableNames', node_index_table_names);

                % beam nodes
                tmp = beam_truss.Kn_ele( [P(S(b1,2),2:4),P(S(b1,3),2:4), n]); 
                for i_b = 1:n-1
                    K.Balken{b1}.el_type{4}.num{i_b}.coordinates = tmp(i_b,:);
                    node_index.beam{b1} = [ node_index.beam{b1}; ...
                        array2table( [b1 , 4 , i_b], 'Variablenames', node_index_table_names ) ];
                end

                % main nodes (end)
                K.Balken{b1}.el_type{1}.num{2}.coordinates = P(S(b1,3),2:4);
                node_index.beam{b1} = [ node_index.beam{b1}; ...
                    array2table([b1 , 1 , 2], 'Variablenames', node_index_table_names) ];

                % metadata table
                node_index.beam{b1}.Properties.VariableDescriptions{'beam'} = ...
                    'number of beam (see matrix Beam)';
                node_index.beam{b1}.Properties.VariableDescriptions{'element_type'} = ...
                    '1 - main node, 4 - beam';
                node_index.beam{b1}.Properties.VariableDescriptions{'node_number'} = ...
                    'consecutive number of the element of the respective element type';
            end

            clear b1

        end

        function K = Kn_ele(in)
            % support function for truss discretization
            % divide the beam into elements and get the coordinates

            xa= in(1);
            ya= in(2);
            za= in(3);
            xe= in(4);
            ye= in(5);
            ze= in(6);
            n = in(7);

            K = zeros(n-1,3);

            for ii = 1:n-1
                K(ii,1) = (xa + (ii/(n-ii))*xe)/(1+(ii/(n-ii)));
                K(ii,2) = (ya + (ii/(n-ii))*ye)/(1+(ii/(n-ii)));
                K(ii,3) = (za + (ii/(n-ii))*ze)/(1+(ii/(n-ii)));
            end
        end

        function [intersectVectors, ind_a, ind_b] = intersectrows_fast(a,b,~)
            % credits: https://stackoverflow.com/questions/28493795/faster-alternative-to-intersect-with-rows-matlab
            % calculate equivalent one-column versions of input arrays

            mult = (10^ceil(log10( 1+max( [a(:);b(:)] ))).^(size(a,2)-1:-1:0))';
            acol1 = a*mult;
            bcol1 = b*mult;

            ind_a=find(acol1==bcol1);
            intersectVectors=[];
            ind_b=[];

        end
       
        function M_glob = M_Matrix(kn_a,kn_e,Ms_loc,C2,B1)
            % subsitute the symbolic variables in the mass and stiffness matrix
            
            syms lsy
            % element length
            l = abs(norm(kn_e-kn_a,2));
            
            % substitute the element length into the matrix
            M_loc=eval(subs(Ms_loc,lsy,l));
            
            % transformation into structural axes
            M_glob=B1.'*C2.'*M_loc* C2*B1;
            
        end    
        
    end

end