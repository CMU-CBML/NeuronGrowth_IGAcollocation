
% IGA-collocation Implementation for 2D neuron growth
% Kuanren Qian
% 10/27/2021

%% CleanUp
close all;
clear;
clc;

% Including Path 
addpath('../IGA_collocation_algorithm');

%% Start Simulation Model
disp('********************************************************************');
disp('2D Phase-field Neuron Growth solver using IGA-Collocation');
disp('********************************************************************');

% save rng seed for repeatability
rngSeed = rng('shuffle');
save('./data/rngSeed','rngSeed');

% variable and png save frequency
var_save_invl = 500;
png_save_invl = 100;
png_plot_invl = 50;

%% Simulation Parameter Initialization
% time stepping variables
dtime = 1.5e-2;
end_iter = 20000;

% setup domain size outside main for modularity (main.m can remain the same
% for different simulation cases)
neuron_domain_setup
 
% B-spline curve order (U,V direction)
p = 3;
q = 3;
% knot vector
knotvectorU = [0,0,0,linspace(0,Nx,Nx+1),Nx,Nx,Nx].';
knotvectorV = [0,0,0,linspace(0,Ny,Ny+1),Ny,Ny,Ny].';

% setting lenu lenv this way for easier access to ghost nodes later on
lenu = length(knotvectorU)-2*(p-1);
lenv = length(knotvectorV)-2*(p-1);

% neuron growth variables
aniso = 6;
kappa= 4;
alph = 0.9; % changing name to alph cause alpha is a function
pix=4.0*atan(1.0);
alphOverPix = alph/pix;
% gamma = 10.0;
gamma = 15.0;
tau = 0.3;
M_phi = 60;
M_theta = 0.5*M_phi;
s_coeff = 0.007;
delta = 0.1;
epsilonb = 0.04;

% Tubulin parameters
r = 1e6;
g = 0.1;
alpha_t = 0.001;
beta_t = 0.001;
Diff = 4;
% source_coeff = 0.02;
source_coeff = 0.05;

% tolerance for NR iterations in phi equation
tol = 1e-4;

% Seed size
seed_radius = 20;
% initializing phi and concentration based on neuron seed position
[phi,conct] = initialize_neurite_growth(seed_radius,lenu,lenv,numNeuron); 

% Expanding domain parameters
BC_clearance = 15;
expd_sz = 10;

disp(' Simulation Parameter Initialization - Done!');

%% Iterating Variable Initialization
% constructing collocation basis
order_deriv = 2;    % highest order of derivatives to calculate
[cm,size_collpts] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv);
lap = cm.N2uNv + cm.NuN2v;
[lap_flip, lap_id] = extract_diags(lap);

phi = cm.NuNv\phi;
conc_t = cm.NuNv\conct;

% initializing theta and temperature
theta=cm.NuNv\reshape(rand(lenu,lenv),lenu*lenv,1);
theta_ori = zeros(lenu,lenv);
tempr = zeros([lenu*lenv,1]);

% theta does not evolve over time, only need to compute initially or expanding domain
% magnitude of theta gradient
mag_grad_theta = sqrt(bsxfun(@times,cm.N1uNv*theta,cm.N1uNv*theta)+bsxfun(@times,cm.NuN1v*theta,cm.NuN1v*theta));
C0 = 0.5+6*s_coeff*mag_grad_theta;

% initializing initial phi,theta,tempr for boundary condition (Dirichlet)
phi_initial = reshape(phi,lenu,lenv);
theta_initial = reshape(theta,lenu,lenv);
tempr_initial = reshape(tempr,lenu,lenv);
for i = 2:lenu-1
    for j = 2:lenv-1
        phi_initial(i,j) = 0;
        theta_initial(i,j) = 0;
        tempr_initial(i,j) = 0;
    end
end
phi_initial = reshape(phi_initial,lenu*lenv,1);
theta_initial  = reshape(theta_initial,lenu*lenv,1);
tempr_initial  = reshape(tempr_initial,lenu*lenv,1);
save('./data/phi_on_cp_initial','phi');
save('./data/theta_on_cp_initial','theta');
save('./data/tempr_on_cp_initial','tempr');

% plotting initial phi
set(gcf,'position',[100,100,800,400]);
colormap parula;

% binary ID for boundary location (define 4 edges)
% id = 1 means there is bc
bcid = zeros([lenu,lenv]);
for i = 1:lenu
    bcid(1,i) = 1;
    bcid(lenu,i) = 1;
    bcid(i,1) = 1;
    bcid(i,lenv) = 1;
end
bcid = reshape(bcid,lenu*lenv,1);

disp('Iterating Variable Initialization - Done!');
disp('********************************************************************');

%% Neuron Growth Stage Variable Initialization
iter_stage2_begin = 500;
iter_stage3_begin = iter_stage2_begin+2000;
iter_stage45_begin = iter_stage3_begin+10000;
rot_iter_invl = 500;

Rot = zeros(1,20);
rotate = zeros(1,20);
rotate_intv = zeros(1,20);

disp('Neuron Growth Stage Variable Initialization - Done!');
disp('********************************************************************');

%% Transient iteration computation
disp('Starting Neuron Growth Model transient iterations...');
tic;
iter = 0;
while iter <= end_iter
    iter = iter + 1;
    if(mod(iter,50) == 0)
        fprintf('Progress: %.2d/%.2d\n',iter,end_iter);
        toc
    end
    
    % calculating a and a*a' (aap) in the equation using theta and phi
    xtheta = cm.NuNv*theta;
    [a, ~, aap,~,~] = kqGetEpsilonAndAap(epsilonb,delta,phi,xtheta,cm.L_NuNv,...
        cm.U_NuNv,cm.NuN1v,cm.N1uNv);
    
    NNtempr = cm.NuNv*tempr;
    NNct = cm.NuNv*conc_t;

    if(iter<=iter_stage2_begin)
        E = alphOverPix*atan(gamma*(1-NNtempr));
    else
        nnT = reshape(theta_ori,lenu*lenv,1);
        delta_L = r*NNct - g;
        term_change = (regular_Heiviside_fun(delta_L));
        term_change(nnT==1)=1;

        E = alphOverPix*atan(gamma*bsxfun(@times,term_change,1-NNtempr));
        E(abs(nnT)==0) = 0;
        
        if(mod(iter,png_plot_invl) == 0)
            subplot(2,3,4);
            phi_plot = reshape(cm.NuNv*phi,lenu,lenv);
            imagesc(reshape(E,lenu,lenv)+phi_plot);
            title(sprintf('E overlay with phi'));
            axis square;
            colorbar;
        end
    end

    %% Phi (Implicit Nonlinear NR method)
    % NR method initial guess (guess current phi)
    phiK = phi;
    % residual for NR method (make sure it is larger than tolerance at the beginning)
    R = 2*tol;

    % splitted C0 from C1 because E mag_grad_theta dimension mismatch
    % during domain expansion. Compute here to fix conflict
    C1 = E-C0;

    % NR method calculation
    ind_check = 0;
    NNa = cm.NuNv*a;
    N1Na = cm.N1uNv*a;
    NN1a = cm.NuN1v*a;
    NNaap = cm.NuNv*aap;
    N1Naap = cm.N1uNv*aap;
    NN1aap = cm.NuN1v*aap;

    out = cell(4,1);
    dt_t = 0;
    while max(abs(R)) >= tol
        NNpk = cm.NuNv*phiK;
        N1Npk = cm.N1uNv*phiK;
        NN1pk = cm.NuN1v*phiK;
        N1N1pk = cm.N1uN1v*phiK;
        LAPpk = lap*phiK;
        
        % term a2
        out{1}  = 2*bsxfun(@times,bsxfun(@times,NNa,N1Na),N1Npk)+bsxfun(@times,bsxfun(@times,NNa,NNa),LAPpk) ...
            +2*bsxfun(@times,bsxfun(@times,NNa,NN1a),NN1pk);
        % termadx
        out{2} = bsxfun(@times,N1Naap,NN1pk)+bsxfun(@times,NNaap,N1N1pk);
        %termady
        out{3} = bsxfun(@times,NN1aap,N1Npk)+bsxfun(@times,NNaap,N1N1pk);
        % termNL
        out{4} = -bsxfun(@times,bsxfun(@times,NNpk,NNpk),NNpk)+bsxfun(@times,bsxfun(@times,(1-C1),NNpk),NNpk)+bsxfun(@times,C1,NNpk);
        if dt_t==0 % these terms only needs to be calculated once
            % terma2_deriv
            t5 =  banded_dot_star((2*bsxfun(@times,NNa,N1Na)+N1Naap),cm.N1uNv_flip,cm.N1uNv_id) +...
                  banded_dot_star(bsxfun(@times,NNa,NNa),lap_flip,lap_id) +...
                  banded_dot_star((2*bsxfun(@times,NNa,NN1a)-N1Naap),cm.NuN1v_flip,cm.NuN1v_id);
        end
        % termNL_deriv
        temp = (-3*bsxfun(@times,NNpk,NNpk)+2*bsxfun(@times,(1-C1),NNpk)+C1);
        t6 = banded_dot_star(temp, cm.NuNv_flip, cm.NuNv_id);
        
        R = M_phi/tau*(out{1}-out{2}+out{3}+out{4});
        R = R*dtime-NNpk+(cm.NuNv*phi);
        dR = M_phi/tau*(t5+t6);
        dR = dR*dtime-cm.NuNv;

        % check residual and update guess
        R = R - dR*phi_initial;
        [dR, R] = StiffMatSetupBCID(dR, R,bcid,phi_initial);
        dp = dR\(-R);
        phiK = phiK + dp;
        
        max_phi_R = full(max(abs(R)));
%         fprintf('Phi NR Iter: %.2d -> max residual: %.2d\n',ind_check, max_phi_R);
        if (ind_check >= 100 || max(abs(R))>1e20)
            error('Phi NR method NOT converging!-Max residual: %.2d\n',max_phi_R);
        end
        ind_check = ind_check + 1;
        dt_t = dt_t+dtime;
    end

    %% Temperature (Implicit method)
    temprLHS = cm.NuNv-3*dt_t*lap;
    temprRHS = kappa*(cm.NuNv*phiK-cm.NuNv*phi)+NNtempr;
    [temprLHS, temprRHS] = StiffMatSetupBCID(temprLHS, temprRHS,bcid,tempr_initial);
    tempr_new = temprLHS\temprRHS;

    %% Tubulin concentration (Implicit method)
    NNp = cm.NuNv*phi;
    nnpk = round(NNpk);

    if iter == 1 % save initial lap phi, this var will not change and will be used throughout the simulation
        initial_LAPpk = lap*phi;
        sum_lap_phi = sum(bsxfun(@times,initial_LAPpk,initial_LAPpk));
        save(sprintf('./data/initial_LAPpk_%2d',iter),'initial_LAPpk');
    end
    
    term_diff = Diff*(banded_dot_star(N1Npk,cm.N1uNv_flip,cm.N1uNv_id) + ...
        banded_dot_star(NNpk,lap_flip,lap_id) + ...
        banded_dot_star(NN1pk,cm.NuN1v_flip,cm.NuN1v_id));
    term_alph = alpha_t*(banded_dot_star(N1Npk,cm.NuNv_flip,cm.NuNv_id) + ...
        banded_dot_star(NNpk,cm.N1uNv_flip,cm.N1uNv_id) + ...
        banded_dot_star(NN1pk,cm.NuNv_flip,cm.NuNv_id) + ...
        banded_dot_star(NNpk,cm.NuN1v_flip,cm.NuN1v_id));
    term_beta = beta_t*banded_dot_star(NNpk,cm.NuNv_flip,cm.NuNv_id);
    term_source = source_coeff/sum_lap_phi*bsxfun(@times,initial_LAPpk,initial_LAPpk);

    conc_t_RHS = dtime/2*term_source-bsxfun(@times,NNct,(NNpk-NNp))+bsxfun(@times,NNp,NNct);
    conc_t_LHS = bsxfun(@times,NNp,cm.NuNv)-dtime/2*(term_diff-term_alph-term_beta);
    bcid_t = (~nnpk);
    [conc_t_LHS, conc_t_RHS] = StiffMatSetupBCID(conc_t_LHS, conc_t_RHS,bcid_t,zeros(lenu*lenv,1));
    conc_t_new = conc_t_LHS\conc_t_RHS;
    
    %% iteration update
    % update variables in this iteration
    phi = phiK;
    tempr = tempr_new;
    conc_t = conc_t_new;

    %% Plotting figures and check for domain expansion
    if(mod(iter,png_plot_invl) == 0 || iter == 1)
        phi_plot = reshape(cm.NuNv*phiK,lenu,lenv);
        tempr_plot = reshape(cm.NuNv*tempr_new,lenu,lenv);
        conct_plot = reshape(cm.NuNv*conc_t_new,lenu,lenv);
        
        subplot(2,3,1);
        imagesc(phi_plot);
        title(sprintf('Phi at iteration = %.2d',iter));
        axis square;
        colorbar;

        subplot(2,3,2);
        imagesc(tempr_plot);
        title(sprintf('T at iteration = %.2d',iter));
        axis square;
        colorbar;

        subplot(2,3,3);
        imagesc(conct_plot);
        title(sprintf('Tubulin at iteration = %.2d',iter));
        axis square;
        colorbar;
        
        % plot current iteration
        drawnow;

        % save picture
        if(mod(iter,png_save_invl) == 0)
            try
                saveas(gcf,sprintf('./data/NeuronGrowth_%.2d.png',iter));
            catch
                fprintf('png write error skipped.\n');
            end
        end
        
        expd_state = 0;
        % check for domain expansion
        if(iter~=1 && (max(max(phi_plot(1:BC_clearance,:))) > 0.5 || ...
                max(max(phi_plot(:,1:BC_clearance))) > 0.5 || ...
                max(max(phi_plot(end-BC_clearance:end,:))) > 0.5 || ...
                max(max(phi_plot(:,end-BC_clearance:end))) > 0.5))
           
            disp('********************************************************************');
            disp('Expanding Domain...');

            Nx = Nx+expd_sz;
            Ny = Ny+expd_sz;

            knotvectorU = [0,0,0,linspace(0,Nx,Nx+1),Nx,Nx,Nx].';
            knotvectorV = [0,0,0,linspace(0,Ny,Ny+1),Ny,Ny,Ny].';

            lenu = length(knotvectorU)-2*(p-1);
            lenv = length(knotvectorV)-2*(p-1);

            oldNuNv = cm.NuNv;
            [cm, size_collpts] = kqCollocationDers(knotvectorU,p,knotvectorV,...
                q,order_deriv);
            lap = cm.N2uNv + cm.NuN2v;
            [lap_flip, lap_id] = extract_diags(lap);

            sz = length(lap);

            [phi,theta,conc_t,tempr,initial_LAPpk,phi_initial,theta_initial,tempr_initial,bcid] ...
                = kqExpandDomain(sz,phiK,theta,conc_t_new,tempr_new,initial_LAPpk,oldNuNv,cm.NuNv);

            phiK = phi;

            mag_grad_theta = sqrt(bsxfun(@times,cm.N1uNv*theta,cm.N1uNv*theta)+bsxfun(@times,cm.NuN1v*theta,cm.NuN1v*theta));
            C0 = 0.5+6*s_coeff*mag_grad_theta;

            expd_state = 1;

            toc
            disp('********************************************************************');
        end
    end
    
    % save variables
    if(mod(iter,var_save_invl)==0 || iter == 0)
        try
            save(sprintf('./data/phi_on_cp_%2d',iter),'phi');
            save(sprintf('./data/tempr_on_cp_%2d',iter),'tempr');
            save(sprintf('./data/conct_on_cp_%2d',iter),'conc_t');
        catch
            fprintf('data write error skipped.\n');
        end
    end
    
    %% Growth stage operations
    if(mod(iter,1000)==0||expd_state==1)
        if iter <=10000
            CX = [lenu,0];
            CY = [0,lenv];
        else
            CX = [lenu,0,0];
            CY = [0,lenv,lenv*2/3];
        end
    end

    phi_plot = full(reshape(cm.NuNv*phiK,lenu,lenv));
    
    if iter>=500 && iter<=1000
        theta_ori = growthStages2(phi_plot);
    elseif iter >= 1000
        [lenu,lenv] = size(phi_plot);
        phi_id = round(phi_plot);
        tip = sum_filter(phi_id,1);
        regionalMaxima = imregionalmax(full(tip));
        L = bwconncomp(regionalMaxima,4);
        S = regionprops(L,'Centroid');

        centroids = floor(cat(1,S.Centroid));
        %     % Matlab MEX cannot directly access S.centroid, need to loop
        %     numRegions = numel(S);
        %     centroids = zeros(numRegions, numel(S(1).Centroid)); % pre-initialize to set size
        %     for region = 1:numRegions
        %         centroids(region, :) = floor(S(region).Centroid); % write one row at a time
        %     end

        max_x = [];
        max_y = [];
        numCue = length(CX);
        cue= zeros(lenu,lenv);
        for k = 1:numCue
            for i = 1:lenu
                for j = 1:lenv
                    cue(i,j) = 1/sqrt(bsxfun(@times,(i-CX(k)),(i-CX(k)))+bsxfun(@times,(j-CY(k)),(j-CY(k))));
                end
            end

            tips = zeros(1,L.NumObjects);
            for i = 1:L.NumObjects
                tips(i) = cue(centroids(i,1),centroids(i,2));
            end
            [~,tip_want_ind] = max(tips);
        %         tip_want_ind = find(tips>=0.95*max(tips));
            max_x = [max_x,centroids(tip_want_ind,1).'];
            max_y = [max_y,centroids(tip_want_ind,2).'];
        end

        size_Max = length(max_x);
        [theta_ori] = theta_rotate(lenu,lenv,max_x,max_y,size_Max);
        %         theta_ori = cuedGrowthStages3(phi_plot,CX,CY);
        %         theta_ori = cuedGrowthStages3_mex(phi_plot,CX,CY);
    end

end

disp('All simulations complete!\n');
