function [T, Q, Qp, niter,Qpp] = kinematic_analysis(mbs, q0, h, t_end, tol)

T = 0.0:h:t_end; %Time (time step h)
nt = length(T);
Q = zeros(nt, mbs.nq);
Qp = zeros(nt, mbs.nq); 
Qpp= zeros(nt,mbs.nq);

qi = q0;
niter = zeros(1, nt); %num of iterations

for ii = 1 : nt
    t = T(ii);
    [qi, niter(ii)] = ...
        NR_method(@(y)constraints(mbs, y, t), ...
        @(y)constraints_dq(mbs, y), qi, tol);
    Q(ii, :) = qi';
    % Below is velocity analysis
    Cq = constraints_dq(mbs, qi); %Cq_fun
    Ct = constraints_dt(mbs, t); %Ct_fun
    qip = -Cq\Ct; % -Cq^-1*Ct
    Qp(ii, :) = qip';
    qi = qi + h .* qip; %reduces number of iterations from 3 to 2
        %old position + time step*velocity
        
        
    %Below is acceleration analysis
    g = constraints_g(mbs,qi,qip,t); 
    qipp = -Cq\g; 
    Qpp(ii,:) = qipp'; 
        
end

