function [phi,conct] = initialize_neurite_growth(seed_radius,lenu,lenv,numNeuron)
% this function is was initially written for sparse matrix. The main.m code
% was later changed to full mat.
seed = (seed_radius)^2;

ind_i = [];
ind_j = [];
phi_val = [];
conct_val = [];

for i=1:lenu
    for j=1:lenv
        if ((i-lenu/2)*(i-lenu/2)+(j-lenv/2)*(j-lenv/2) < seed)
            r = sqrt((i-lenu/2)*(i-lenu/2)+(j-lenv/2)*(j-lenv/2));
            if numNeuron == 1 % 1 neuron
                ind_i(end+1) = i;
                ind_j(end+1) = j;
                phi_val(end+1) = 1;            
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
            elseif numNeuron == 2 % 2 neurons
                ind_i(end+1) = i+90;
                ind_i(end+1) = i-90;
                ind_j(end+1) = j+90;
                ind_j(end+1) = j-90;
                
                phi_val(end+1) = 1;            
                phi_val(end+1) = 1;                    
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
            elseif numNeuron == 3 % 3 neurons
                ind_i(end+1) = i-60;
                ind_i(end+1) = i+60;
                ind_i(end+1) = i+60;
                
                ind_j(end+1) = j;
                ind_j(end+1) = j+40;
                ind_j(end+1) = j-40;  
                
                phi_val(end+1) = 1;            
                phi_val(end+1) = 1;                    
                phi_val(end+1) = 1;      
                
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
            elseif numNeuron == 4 % 4 neurons
                ind_i(end+1) = i-60;
                ind_i(end+1) = i+60;
                ind_i(end+1) = i+60;
                ind_i(end+1) = i-60;
                
                ind_j(end+1) = j-60;
                ind_j(end+1) = j+60;
                ind_j(end+1) = j-60;  
                ind_j(end+1) = j+60;  
                
                phi_val(end+1) = 1;            
                phi_val(end+1) = 1;                    
                phi_val(end+1) = 1;      
                phi_val(end+1) = 1;      
                
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
            elseif numNeuron == 5 % 5 neurons
                ind_i(end+1) = i-75;
                ind_i(end+1) = i+75;
                ind_i(end+1) = i+75;
                ind_i(end+1) = i-75;
                ind_i(end+1) = i;
                
                ind_j(end+1) = j-75;
                ind_j(end+1) = j+75;
                ind_j(end+1) = j-75;  
                ind_j(end+1) = j+75;  
                ind_j(end+1) = j;  
                
                phi_val(end+1) = 1;            
                phi_val(end+1) = 1;                    
                phi_val(end+1) = 1;      
                phi_val(end+1) = 1;      
                phi_val(end+1) = 1;      
                
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
            elseif numNeuron == 6 % 6 neurons
                ind_i(end+1) = i-85;
                ind_i(end+1) = i-10;
                ind_i(end+1) = i+85;
                ind_i(end+1) = i-85;
                ind_i(end+1) = i+10;
                ind_i(end+1) = i+85;
                
                ind_j(end+1) = j-15;
                ind_j(end+1) = j-85;
                ind_j(end+1) = j-85;  
                ind_j(end+1) = j+85;  
                ind_j(end+1) = j+85;  
                ind_j(end+1) = j+15;  
                
                phi_val(end+1) = 1;            
                phi_val(end+1) = 1;                    
                phi_val(end+1) = 1;      
                phi_val(end+1) = 1;      
                phi_val(end+1) = 1;      
                phi_val(end+1) = 1;      
                
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
            elseif numNeuron == 7 % 7 neurons
                ind_i(end+1) = i-95;
                ind_i(end+1) = i-45;
                ind_i(end+1) = i+45;
                ind_i(end+1) = i-45;
                ind_i(end+1) = i+45;
                ind_i(end+1) = i+95;
                ind_i(end+1) = i;

                ind_j(end+1) = j;
                ind_j(end+1) = j-95;
                ind_j(end+1) = j-95;  
                ind_j(end+1) = j+95;  
                ind_j(end+1) = j+95;  
                ind_j(end+1) = j;  
                ind_j(end+1) = j;  

                phi_val(end+1) = 1;            
                phi_val(end+1) = 1;                    
                phi_val(end+1) = 1;      
                phi_val(end+1) = 1;      
                phi_val(end+1) = 1;      
                phi_val(end+1) = 1;      
                phi_val(end+1) = 1;      
                
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
                conct_val(end+1) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
            
            end
        end
    end
end

%% Creating sparse matrix
phi = sparse(ind_i,ind_j,phi_val,lenu,lenv);
phi = full(reshape(phi,lenu*lenv,1)); % reshpae phi for calculation
conct = sparse(ind_i,ind_j,conct_val,lenu,lenv);
conct = full(reshape(phi,lenu*lenv,1)); % reshpae phi for calculation
