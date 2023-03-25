clear
clc
close all
% T,Torque and Force are introduced from Book1.xlsx,
% saved in Yin.mat
load('input_data.mat');
% By using ginput function,
% each start and end point of the drilling process 
% is saved in yinkexin.mat
load('Split_points.mat');
% Set sample Frequency
Fs = 10000; 
%% Calibration
% Since the raw signals we get from LabView are counts,
% methods should be used to 
% convert counts to true forces and torques via voltage. 
% Convert count to voltage
VOF = (20/4096).*Force; 
VOT = (20/4096).*Torque;
% Convert voltage to force and torque
TrueForce = 1000.*VOF;
TrueTorque = 10.*VOT;
TrueT=t/10^6;
%% Raw Signals
for i=1:9
SingleTorque=TrueTorque(yinkexin(2*i-1):yinkexin(2*i));
SingleTime=TrueT(yinkexin(2*i-1):yinkexin(2*i));
yyaxis left;
ylabel('Torque (Nm)')
ylim([-10,40]);
subplot(5,2,i);
plot(SingleTime,SingleTorque);
hold on;
SingleForce=TrueForce(yinkexin(2*i-1):yinkexin(2*i))/1000;
SingleTime=TrueT(yinkexin(2*i-1):yinkexin(2*i));
xlabel('Time (s)')
yyaxis right;
ylabel('Force (kN)')
ylim([0,10]);
plot(SingleTime,SingleForce,'r');
grid on;
end
%% Decimation
Block=1000;
for j=1:9
BlockNumber=floor(length(TrueTorque(yinkexin(2*j-1):yinkexin(2*j)))/Block);
BlockTorque=ones(1,BlockNumber);
BlockForce=ones(1,BlockNumber);
BlockT=ones(1,BlockNumber);
BlockStart=floor((length(0:yinkexin(2*j-1)))/Block);
BlockEnd=floor((length(0:yinkexin(2*j)))/Block);
  for i=BlockStart:BlockEnd
    BlockTorque(i)=median(TrueTorque(Block*(i-1)+1:Block*i));
    subplot(5,2,j);
    plot(BlockTorque,'- .','MarkerSize',6,'LineWidth',1);
    hold off
    hold on
    yyaxis left
    ylabel('Torque (Nm)');
    xlabel('Number of Blocks');
    grid on;
    xlim([BlockStart,BlockEnd]);
    BlockForce(i)=median((TrueForce(Block*(i-1)+1:Block*i)))/1000;
    plot(BlockForce,'- .','MarkerSize',6,'LineWidth',1);
    yyaxis right;
    ylabel('Force (kN)');
    hold off
  end
end
%% Remove outliers
for i=1:9
figure;
set_1=TrueTorque(yinkexin(2*i-1):yinkexin(2*i));
T=(t(yinkexin(2*i-1):yinkexin(2*i)))/10^6;
plot(T,set_1)
hold on;
star_1=filloutliers(set_1,'nearest','mean');
plot(T,star_1);
xlabel('Time (s)')
ylabel('Torque (Nm)')
legend('Raw data','Remove Outliers')
grid on;
end
for i=1:9
figure;
set_2=TrueForce(yinkexin(2*i-1):yinkexin(2*i));
T=(t(yinkexin(2*i-1):yinkexin(2*i)))/10^6;
plot(T,set_2)
hold on;
star_2=filloutliers(set_2,'nearest','mean');
plot(T,star_2);
xlabel('Time (s)')
ylabel('Torque (Nm)')
end
%% Max Torque
MaxTorque=ones(1,9);
for i=1:9
set_1=TrueTorque(yinkexin(2*i-1):yinkexin(2*i));
star_1=filloutliers(set_1,'nearest','mean');
MaxTorque(i)=max(star_1);
end
%the effects of feed rate on max torque
Group_MaxTorque_1=[MaxTorque(1),MaxTorque(2),MaxTorque(3);MaxTorque(4),MaxTorque(5),MaxTorque(6);MaxTorque(7),MaxTorque(8),MaxTorque(9)];
bar(Group_MaxFTorque_1,1,'grouped')
ylim([0,30]);
ylabel('Max Torque (Nm)')
xlabel('Spindle speed (rpm)')
set(gca,'xticklabel',1114:477:2068)
legend('96 mm/min','137 mm/min','178 mm/min')
grid on
%the effects of speed on max torque
Group_MaxTorque_2=[MaxTorque(1),MaxTorque(4),MaxTorque(7);MaxTorque(2),MaxTorque(5),MaxTorque(8);MaxTorque(3),MaxTorque(6),MaxTorque(9)];
bar(Group_MaxTorque_2,1,'grouped')
ylim([0,32])
ylabel('Max Torque (Nm)')
xlabel('Feed Rate (mm/min)')
set(gca,'xticklabel',96:41:178)
legend('1114 rpm','1591 rpm','2068 rpm')
grid on
%% Mean Torque
MeanTorque=ones(1,9);
for i=1:9
set_1=TrueTorque(yinkexin2(2*i-1):yinkexin2(2*i));
star_1=filloutliers(set_1,'nearest','mean');
MeanTorque(i)=mean(star_1);
end
%the effects of feed rate on mean torque
Group_MeanTorque_1=[MeanTorque(1),MeanTorque(2),MeanTorque(3);MeanTorque(4),MeanTorque(5),MeanTorque(6);MeanTorque(7),MeanTorque(8),MeanTorque(9)];
bar(Group_MeanTorque_1,1,'grouped')
ylim([0,14]);
ylabel('Mean Torque (Nm)')
xlabel('Spindle speed (rpm)')
set(gca,'xticklabel',1114:477:2068)
legend('96 mm/min','137 mm/min','178 mm/min')
grid on
%the effects of speed on mean torque
Group_MeanTorque_2=[MeanTorque(1),MeanTorque(4),MeanTorque(7);MeanTorque(2),MeanTorque(5),MeanTorque(8);MeanTorque(3),MeanTorque(6),MeanTorque(9)];
bar(Group_MeanTorque_2,1,'grouped')
ylim([0,16]);
ylabel('Mean Torque (Nm)')
xlabel('Feed Rate (mm/min)')
set(gca,'xticklabel',96:41:178)
legend('1114 rpm','1591 rpm','2068 rpm')
grid on
%% Max Force
MaxForce=ones(1,9);
for i=1:9
set_1=TrueForce(yinkexin(2*i-1):yinkexin(2*i))/1000;
star_1=filloutliers(set_1,'nearest','mean');
MaxForce(i)=max(star_1);
end
%the effects of feed rate on max force
Group_MaxForce_1=[MaxForce(1),MaxForce(2),MaxForce(3);MaxForce(4),MaxForce(5),MaxForce(6);MaxForce(7),MaxForce(8),MaxForce(9)];
bar(Group_MaxForce_1,1,'grouped')
grid on
ylim([0,12])
ylabel('Max Force (kN)')
xlabel('Spindle speed (rpm)')
set(gca,'xticklabel',1114:477:2068)
legend('96 mm/min','137 mm/min','178 mm/min')
grid on
%the effects of speed on max force
Group_MaxForce_2=[MaxForce(1),MaxForce(4),MaxForce(7);MaxForce(2),MaxForce(5),MaxForce(8);MaxForce(3),MaxForce(6),MaxForce(9)];
bar(Group_MaxForce_2,1,'grouped')
ylim([0,13]);
ylabel('Max Force (kN)')
xlabel('Feed Rate (mm/min)')
set(gca,'xticklabel',96:41:178)
legend('1114 rpm','1591 rpm','2068 rpm')
grid on
%% Mean Force
MeanForce=ones(1,9);
for i=1:9
set_1=(TrueForce(yinkexin(2*i-1):yinkexin(2*i))/1000);
star_1=filloutliers(set_1,'nearest','mean');
MeanForce(i)=mean(star_1);
end
%the effects of feed rate on mean force
Group_MeanForce_1=[MeanForce(1),MeanForce(2),MeanForce(3);MeanForce(4),MeanForce(5),MeanForce(6);MeanForce(7),MeanForce(8),MeanForce(9)];
bar(Group_MeanForce_1,1,'grouped')
grid on
ylim([0,8])
ylabel('Mean Force (kN)')
xlabel('Spindle speed (rpm)')
set(gca,'xticklabel',1114:477:2068)
legend('96 mm/min','137 mm/min','178 mm/min')
grid on
%the effects of speed on mean force
Group_MeanForce_2=[MeanForce(1),MeanForce(4),MeanForce(7);MeanForce(2),MeanForce(5),MeanForce(8);MeanForce(3),MeanForce(6),MeanForce(9)];
bar(Group_MeanForce_2,1,'grouped')
ylim([0,9])
ylabel('Mean Force (kN)')
xlabel('Feed Rate (mm/min)')
set(gca,'xticklabel',96:41:178)
legend('1114 rpm','1591 rpm','2068 rpm')
grid on
%% Fitting
fym=[0.1407	0.1870	0.2305	0.1058	0.1406	0.1733	0.0858	0.1140	0.1406];
fyf=[0.1798	0.2306	0.2770	0.1400	0.1797	0.2158	0.1166	0.1496	0.1796];
%% Detrending the signal 
% citing from Vallabha Hampiholi,
% using function movavgFilt allows us to detrend the signal 
% calculate the difference of true torque and moving average.
i=1;
ATorque = (TrueTorque(yinkexin(2*i-1):yinkexin(2*i)))';
T1=(t(yinkexin(2*i-1):yinkexin(2*i)))/10^6;
movavgTorque = movavgFilt(ATorque,101,'Center');
ATorque = ATorque-movavgTorque;
figure
subplot(2,1,1)
hold on
plot(T1,movavgTorque)
plot(T1,ATorque,'r-')
hold off 
xlabel('Time (s)')
ylabel('Torque(Nm)')
% calculate the difference of true force and moving average. 
AForce = (TrueForce(yinkexin(2*i-1):yinkexin(2*i)))';
T2=(t(yinkexin(2*i-1):yinkexin(2*i)))/10^6;
movavgTorque = movavgFilt(AForce,101,'Center');
AForce = AForce-movavgTorque;
subplot(2,1,2)
hold on
plot(T2,movavgTorque)
plot(T2,AForce,'r-')
hold off
xlabel('Time (s)')
ylabel('Force(N)')
%% Fast Fourier Transform
% apply FFT on the Torque
A = abs(fft(ATorque)/length(ATorque));
Amp = A(1:length(ATorque)/2+1);
Amp(2:end-1) = 2*Amp(2:end-1);
f = Fs*(0:(length(ATorque)/2))/length(ATorque);
figure
plot(f,Amp)
xlabel('f (Hz)')
ylabel('|Amp(f)|')
grid on;
% apply FFT on the Force
A = abs(fft(AForce)/length(AForce));
Amp = A(1:length(AForce)/2+1);
Amp(2:end-1) = 2*Amp(2:end-1);
f = Fs*(0:(length(AForce)/2))/length(AForce);
figure
plot(f,Amp)
xlabel('f (Hz)')
ylabel('|Amp(f)|')
grid on;
%% 
 length(yinkexin(1):yinkexin(2))