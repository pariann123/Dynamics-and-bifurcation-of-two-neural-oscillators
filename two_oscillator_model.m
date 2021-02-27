function two_oscillator_model
%% Parameters
%portion of firing neuronsin eachpopulation
E1 = 0.5;
I1 = 0.5;
E2 = 0.5;
I2 = 0.5;
% connection strenghts between excitatory and inhibitory populations of SAME
% subpopulation
c1 = 16;
c2 = 12;
c3 = 15;
c4 = 3;
% strengths of the connections btwn DIFFERENT subpopulations related to different
% oscillators
alpha1 = 10;
alpha2 = 0;
alpha3 = 0;
alpha4 = 0;
% external inputs to inhibitory and excitatory populations
Q = 0;
P = 1.5;

theta_e = 4;
theta_i = 3.7;
bi = 2;
be = 1.3;
ke = 1-(1/exp(be*theta_e));
ki = 1-(1/exp(bi*theta_i));

%% Functions
    function dxdt = f(t,x)
        E1 = x(1);
        I1 = x(2);
        E2 = x(3);
        I2 = x(4);
        % connections between the oscillators 
        P12 = alpha1*E2 - alpha2*I2;
        Q12 = alpha3*E2 - alpha4*I2;
        P21 = alpha1*E1 - alpha2*I1;
        Q21 = alpha3*E1 - alpha4*I2;
        dxdt = [-E1+((ke-E1)*(S(c1*E1-c2*I1+P+P12,be,theta_e)));
            -I1+((ki-I1)*(S(c3*E1-c4*I1+Q+Q12,bi,theta_i)));
            -E2+((ke-E2)*(S(c1*E2-c2*I2+P+P21,be,theta_e)));
            -I2+((ki-I2)*(S(c3*E2-c4*I2+Q+Q21,bi,theta_i)))];
    end

    function sigmoid = S(x,b,theta)
        sigmoid = (1./(1+exp(-b*(x-theta))))-(1/(exp(b*theta)));
    end

%% steady state
    function outf= ff(x)
        E = x(1);
        I = x(2);
        % connections between the oscillators 
        P12 = alpha1*E - alpha2*I;
        Q12 = alpha3*E - alpha4*I;
        outf = [-E+((ke-E)*(S(c1*E-c2*I+P+P12,be,theta_e)));
                -I+((ki-I)*(S(c3*E-c4*I+Q+Q12,bi,theta_i)))];
    end


x0s = [0.46,0.49; 0.22,0.15; 0.16,0.12;0.25,0.22;0.02,0.3];
for i=1:length(x0s)
    solve(i,:) = fsolve(@ff, x0s(i,:))
    
end
ss=solve(1,:)
ss2=solve(2,:)


%%  Plot
% initialise figures and clear them
% figure(1)
% clf
% figure(2)
% clf
figure(3)
clf
figure(4)
clf
figure(5)
clf

% set the starting positions
% x0 = [0.3,0.6,0.3,0.6;0.22,0.15,0.25,0.22; 0.16,0.12,0.25,0.22];
% x0 = [0.4,0.5,0.3,0.5]


%%SE, two initial conditions near s.s
% x0 = [0.3,0.6,0.3,0.6];

% x0 = [0.22,0.15,0.16,0.12];
% x0 = [0.16,0.12,0.25,0.2];

%%NC1 and NC2 - one initial condition at S.s the other near limit cycle
x0 = [0.22,0.15,0.25,0.22];
% x0 = [0.2,0.7,0.2,0.7];

%%C and C, both near limit cycle at diff points 
% x0 = [0.16,0.12,0.25,0.22];

% soltuions for E1(t), I1(t), E2(t) and I2(t)
for i = 1:size(x0,1)
    tspan = [0 1000];
    % x gives us 4 columns of each coordinate with rows of ......
    [t,x] = ode45(@f,tspan,x0(i,:,:,:));
    
%     figure(1)
%     subplot(3,1,i)
%     plot(t,x)
%     hold on
%     legend('E1','I1','E2','I2')
%     title('solutions')
%     xlabel('E1(t),I1(t),E2(t),I2(t)')
%     ylabel('Starting positions')
%     
    
    %E1 and E2 on seperate graphs
    figure(5)
    subplot(4,1,1)
    plot(t,x(:,1))
    xlabel('Time, t)')
    ylabel('E1(t)')
    subplot(4,1,2)
    plot(t,x(:,1))
    hold on
    plot(t,x(:,3))
%     xlabel('Time, t','Fontsize',16)
%     ylabel('E1(t),E2(t)','Fontsize',16)
%     legend('E1','E2','Fontsize',12)
%     annotation('textbox',[.9 .5 .1 0.2],'String','\alpha_3=0.1','EdgeColor','none','Fontsize',12)

    subplot(4,1,3)
    plot(t,x(:,1))
    hold on
    plot(t,x(:,4))
    xlabel('Time, t')
    ylabel('E1(t),I2(t)')
    legend('E1','I2')
    subplot(4,1,4)
    plot(t,x(:,3))
    xlabel('Time, t')
    ylabel('E2(t)')



%     phase plane diagram - shows the whole time
%     figure(2)
%     it = floor(length(t)/2);
%     plot(x(it:end,2),x(it:end,4));
%     hold on
%     title('phase-plane diagram')
%     xlabel('I1')
%     ylabel('I2')
    
%     Phase-plane diagrams
    figure(3)
    it = floor(length(t)/2);
    plot(x(it:end,1),x(it:end,3))
    hold on
%     plot(x(it:end,2),x(it:end,4))
    xlabel('E1','fontsize',60)
%     ylabel('E2','fontsize',60)

    figure(4)
    it = floor(length(t)/2);
    it = 1;
    plot(x(it:end,1),x(it:end,2))
    hold on
    plot(x(it:end,3),x(it:end,4))
    xlabel('E','Fontsize',40)
    ylabel('I','Fontsize',40)
%     hold on
%     plot(ss(1,1),ss(1,2),'*','linewidth',5)
    hold on
    plot(ss2(1,1),ss2(1,2),'*','linewidth',5)
end



end
