if (iter >= actin_start)
    len = length(NNpk);
    pp = full(NNpk);
    for i=1:len
        if(pp(i)>0.05)
            pp(i) = 1;
        else
            pp(i) = 0;
        end
    end
    param.RegionMask = reshape(pp,lenu,lenv);
    param.phi = reshape(pp,len,1);

    % 2) set up information for ODE solver
    T   = [0 param.dt];
    sz  = size(param.A);

    Y_0 = [reshape(param.A,sz(1)*sz(2),1); reshape(param.H,sz(1)*sz(2),1)];

    % 3) call ODE solver
    [~,Y] = ode45(@evolvesystem,T,Y_0,[],param);

    sz = size(param.H);
    param.A = reshape(Y(end,1:end/2),sz);
    param.H = reshape(Y(end,end/2+1:end),sz);
    param.T_current = iter*param.dt;

    %% 4) update Hem
    % 1) Change from state 0 to 1
    H_C = max(0,(param.HTotal - sum(sum(param.H)))/param.HTotal);
    kSA = conv2(param.H,param.G_SA,'same');
    lmda = H_C*param.dt*kSA;

    Q = poissrnd(lmda,size(param.H));
    len = sqrt(length(param.phi));
    mk1 = (param.HState == 0) & Q & (param.RegionMask) & reshape(param.phi,len,len);

    % 2) Change from state 1 to 0
    mk0 = (param.HState == 1) & (param.H < param.H_thresh_off) & (param.RegionMask);

    % 3) Update
    % reset param values
    param.H(mk0) = 0;
    param.H(mk1) = param.H_thresh_on;
    param.A(mk0) = 0;

    % reset states
    param.HState(mk0) = 0;
    param.HState(mk1) = 1;
end