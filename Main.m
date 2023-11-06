% 日期：23.7.24
% 作者：Ruixin Miao
ni = 1;
tdd = logspace(-6,-3,3);
t1 = zeros(1,length(tdd));
t2 = zeros(1,length(tdd));
forwarding.MB0.Hsum = 0.0094;
forwarding.model.WT1 = 5;
forwarding.model.WT2 = 1;
forwarding.Pulse.t_90 = 7.9e-5;
forwarding.TLoop.fT = 4e5;
forwarding.TLoop.turn = 1;
forwarding.TLoop.I_trans = 1;
forwarding.B0.Hsum = forwarding.MB0.Hsum;
B1 = 2.5951e-4;
for i = 1:length(tdd)
    i
    forwarding.Pulse.t_90 = tdd(i);
    tic
    [M0_Iter] = IntegralM0_IterM(forwarding,B1,ni);
    t2(i) = toc;
    tic
    [M0_G] = IntegralM0_GM(forwarding,B1,ni,M0_Iter.t);
    t1(i) = toc;
    LGt = length(M0_Iter.t);
    LIt = length(M0_G.t);
    xxx = abs(M0_G.rx - M0_Iter.rx)./abs(M0_Iter.rx);
    xxx(isnan(xxx)) = 0;
    yyy = abs(M0_G.ry - M0_Iter.ry)./abs(M0_Iter.ry);
    yyy(isnan(yyy)) = 0;
    zzz = abs(M0_G.rz - M0_Iter.rz)./abs(M0_Iter.rz);
    zzz(isnan(zzz)) = 0;
    mape_x(i) = mean(xxx,'all');
    mape_y(i) = mean(yyy,'all');
    mape_z(i) = mean(zzz,'all');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% xyz
%%%%%%%%%%%%% x
figure
plot(M0_Iter.t,M0_Iter.rx,'rd:','MarkerIndices',1:5:length(M0_Iter.t),'MarkerSize',1)
hold on
plot(M0_Iter.t,M0_G.rx,'gs-.','MarkerIndices',1:5:length(M0_Iter.t),'MarkerSize',1)
hold off
title('Magnetization vector x')
legend('Geometry analysis','4th-order 5th-order Runge-Kutta')
xlabel('Pulse Emission Time/s');
ylabel('Magnetization/A/m');
%%%%%%%%%%%%% y
figure
plot(M0_Iter.t,M0_Iter.ry,'rd:','MarkerIndices',1:5:length(M0_Iter.t),'MarkerSize',1)
hold on
plot(M0_Iter.t,M0_G.ry,'gs-.','MarkerIndices',1:5:length(M0_Iter.t),'MarkerSize',1)
hold off
title('Magnetization vector y')
legend('Geometry analysis','4th-order 5th-order Runge-Kutta')
xlabel('Pulse Emission Time/s');
ylabel('Magnetization/A/m');
%%%%%%%%%%%%% z 
figure
plot(M0_Iter.t,M0_Iter.rz,'rd:','MarkerIndices',1:5:length(M0_Iter.t),'MarkerSize',1)
hold on
plot(M0_Iter.t,M0_G.rz,'gs-.','MarkerIndices',1:5:length(M0_Iter.t),'MarkerSize',1)
hold off
title('Magnetization vector z')
legend('Geometry analysis','4th-order 5th-order Runge-Kutta')
xlabel('Pulse Emission Time/s');
ylabel('Magnetization/A/m');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAPE
figure
plot(tdd,mape_x,'rd:','MarkerFaceColor','r')
hold on
plot(tdd,mape_y,'gs-.','MarkerFaceColor','g')
plot(tdd,mape_z,'b^--','MarkerFaceColor','b')
hold off
title('MAPE')
legend('MAPE_x','MAPE_y','MAPE_z')
ax5 = gca;
set(ax5,'YScale','log','XScale','log','XGrid','on','YGrid','on')
xlabel('Pulse Emission Time/s');
ylabel('MAPE');
figure
plot(tdd,t1,'rd:','MarkerFaceColor','r')
hold on
plot(tdd,t2,'b^--','MarkerFaceColor','b')
hold off
title('time consuming')
legend('Geometry analysis time-consuming ','4th-order 5th-order Runge-Kutta time-consuming')
ax3 = gca;
set(ax3,'YScale','log','XScale','log','XGrid','on','YGrid','on')
xlabel('Pulse Emission Time/s');
ylabel('Time Consuming/s');