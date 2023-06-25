%% PLANNER (ASSUMING THE POINTS IN path HAVE ALREADY BEEN GENERATED WITH RRT)
Ts = 1e-3;
npoints = length(path(:, 1));
%npoints = 3;
tseg = 5;
ttot=(npoints-1)*tseg; %duration
tdead=10; %dead time to evaluate the steady-state
tot_time = ttot + tdead;

t1=linspace(0,tseg,round(tseg/Ts));
t=linspace(0,ttot+tdead,round(ttot/Ts)+round(tdead/Ts));



%p_d= zeros(4,length(t)); dot_p_d= zeros(4,length(t)); ddot_p_d= zeros(4,length(t));
p_d = []; dot_p_d = []; ddot_p_d = [];

for i = 1:npoints-1
    pi= zeros(4,length(t1)); dot_pi= zeros(4,length(t1)); ddot_pi= zeros(4,length(t1));
    %Initial and final conditions
    x0= path(i, 1); xf=path(i+1, 1); dot_x0= 0; dot_xf=0; ddot_x0=0; ddot_xf=0; dddot_x0 = 0; dddot_xf = 0;
    y0= path(i, 2); yf=path(i+1, 2); dot_y0= 0; dot_yf=0; ddot_y0=0; ddot_yf=0; dddot_y0 = 0; dddot_yf = 0;
    z0= path(i, 3); zf=path(i+1, 3); dot_z0= 0; dot_zf=0; ddot_z0=0; ddot_zf=0; dddot_z0 = 0; dddot_zf = 0;
    psi0 = atan2(y0, x0);
    psif=atan2(yf, xf); dot_psi0= 0; dot_psif=0; ddot_psi0=0; ddot_psif=0; dddot_psi0 = 0; dddot_psif = 0;
    p0=[x0,y0,z0,psi0]; dot_p0=[dot_x0,dot_y0,dot_z0,dot_psi0]; ddot_p0=[ddot_x0,ddot_y0,ddot_z0,ddot_psi0]; dddot_p0=[dddot_x0,dddot_y0,dddot_z0,dddot_psi0];
    pf=[xf,yf,zf,psif]; dot_pf=[dot_xf,dot_yf,dot_zf,dot_psif]; ddot_pf=[ddot_xf,ddot_yf,ddot_zf,ddot_psif]; dddot_pf=[dddot_xf,dddot_yf,dddot_zf,dddot_psif];
    %5-th order polynomial
    a0=zeros(1,4); a1=zeros(1,4); a2=zeros(1,4); a3=zeros(1,4); a4=zeros(1,4); a5=zeros(1,4); a6 = zeros(1,4); a7 = zeros(1,4);
    
    t_iniz = 0;
    t_fin = tseg;

    for j=1:4
        A = [t_iniz^7, t_iniz^6, t_iniz^5, t_iniz^4, t_iniz^3, t_iniz^2, t_iniz, 1;
            t_fin^7, t_fin^6, t_fin^5, t_fin^4, t_fin^3, t_fin^2, t_fin, 1;
            7*t_iniz^6, 6*t_iniz^5, 5*t_iniz^4, 4*t_iniz^3, 3*t_iniz^2, 2*t_iniz, 1, 0;
            7*t_fin^6, 6*t_fin^5, 5*t_fin^4, 4*t_fin^3, 3*t_fin^2, 2*t_fin, 1, 0;
            42*t_iniz^5, 30*t_iniz^4, 20*t_iniz^3, 12*t_iniz^2, 6*t_iniz, 2, 0, 0;
            42*t_fin^5, 30*t_fin^4, 20*t_fin^3, 12*t_fin^2, 6*t_fin, 2, 0, 0;
            210*t_iniz^4, 120*t_iniz^3, 60*t_iniz^2, 24*t_iniz, 6, 0, 0, 0;
            210*t_fin^4, 120*t_fin^3, 60*t_fin^2, 24*t_fin, 6, 0, 0, 0];
        b = [p0(j) pf(j) dot_p0(j) dot_pf(j) ddot_p0(j) ddot_pf(j) dddot_p0(j) dddot_pf(j)]';
        a_temp = A\b;
        a7(j) = a_temp(1);
        a6(j) = a_temp(2);
        a5(j) = a_temp(3);
        a4(j) = a_temp(4);
        a3(j) = a_temp(5);
        a2(j) = a_temp(6);
        a1(j) = a_temp(7);
        a0(j) = a_temp(8);
        
        %trajectories
        pi(j,:)=a7(j)*t1.^7 + a6(j)*t1.^6 + a5(j)*t1.^5 +a4(j)*t1.^4 +a3(j)*t1.^3 +a2(j)*t1.^2 +a1(j)*(t1-(i-1)*tseg) +a0(j);
        dot_pi(j,:) = 7*a7(j)*t1.^6 + 6*a6(j)*t1.^5 + 5*a5(j)*t1.^4 +4*a4(j)*t1.^3 +3*a3(j)*t1.^2 +2*a2(j)*(t1-(i-1)*tseg) +a1(j);
        ddot_pi(j,:) = 42*a7(j)*t1.^5 + 30*a6(j)*t1.^4 + 5*4*a5(j)*t1.^3 +4*3*a4(j)*t1.^2 +3*2*a3(j)*(t1-(i-1)*tseg) +2*a2(j);
    end

    % figure(i)
    % plot(t1, pi)

    p_d = [p_d pi];
    dot_p_d = [dot_p_d dot_pi];
    ddot_p_d = [ddot_p_d ddot_pi];

end

%addition of the stead-state terms
p_d=[p_d zeros(1,round(tdead/Ts))+p_d(:, end)]; 
dot_p_d=[dot_p_d zeros(1,round(tdead/Ts))+dot_p_d(:, end)];
ddot_p_d=[ddot_p_d zeros(1,round(tdead/Ts))+ddot_p_d(:, end)];