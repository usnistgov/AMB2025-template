% Alexander Osterbaan - Bowman Lab - Jan 2025
% Import and plot combined FT-IR (from .srs) and Rheometer (from .xls)
clc;clear;close all

%% Manual inputs
LightOnTime = 40; % (s)
endTime = 55; % (s)

n = 22836; %See FTIR_IndexImport_Seeker.m -- Need this to translate binary from correct starting point

% Peak(s) of interest (will process with linear baseline)
lower1 = 1604;
upper1 = 1651;

lower2 = 1373;
upper2 = 1430;

lower3 = 794;
upper3 = 829;


% --------------- General --------------- 
LineThickness = 3;
FS=12; % font size
figWidth = 540; %540 standard
figHeight = 400; %400 standard
Colors = ["#000000","#e69f00","#56b4e9","#009e73","#f0e442","#0072b2","#d55e00","#cc79a7"]; %Colorblind friendly palette
legendLabel = {'Conversion', 'Light On'};
LegendPos = 'southeast';

% axis variables
xlab = "Time (s)";
xlow = 0;
xlab2 = "Wavenumber (1/cm)";
ylab2 = 'Absorbance (arb. units)';

%======================================================================================================
%% Auto-Inputs
% ID FT-IR
[file,path] = uigetfile('*.srs');
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end

OT = file(12:end-4); %Output title beginning name

% Filename for FT-IR Spectrum Output
OutputTitleFTIR1 = append(OT,'_Peak1_FT-IR.tiff');
OutputTitleFTIR2 = append(OT,'_Peak2_FT-IR.tiff');
OutputTitleFTIR3 = append(OT,'_Peak3_FT-IR.tiff');
% Filename for Kinetics
OutputTitleKinetics1 = append(OT,'_Kinetics1.tiff');
OutputTitleKinetics2 = append(OT,'_Kinetics2.tiff');
OutputTitleKinetics3 = append(OT,'_Kinetics3.tiff');
% Filename for xlsx's
OutputTitle1 = append(OT,'_Kinetics.csv');
OutputTitleRawTime = append(OT,'_RawTime.csv');
OutputTitleRaw= append(OT,'_RawData.csv');

%% Import from FT-IR
% Adapted from version for Omnic 6:  Kurt Oldenburg - 01/11/18, ensure license in folder
fp = fopen(fullfile(path,file));
fseek(fp,14036,'bof');
Spectrum_Pts = fread(fp,1,'int32');
fseek(fp,14048,'bof');
Max_Wavenum=fread(fp,1,'single');
Min_Wavenum=fread(fp,1,'single');
Interval=(Max_Wavenum-Min_Wavenum)/(Spectrum_Pts);
Wavenumbers=linspace(Min_Wavenum,Max_Wavenum,Spectrum_Pts).';
Wavenumbers=flipud(Wavenumbers);
fseek(fp,n,'bof'); % n appears to be dependent on settings, including resolution and number of scans
i=1;
Scratch=fread(fp,25,'uint'); % Spectrum header
Time (i,:)=Scratch(2)/6000;   % Column array of spectra times in minutes
SpectraData (i,:)=fread(fp,Spectrum_Pts,'single');
while ~feof(fp)
    Scratch=fread(fp,25,'uint'); % Spectrum header
    if (Scratch(2)/6000)<Time(i)
    break
    end
    i=i+1; 
    Time (i,:)=Scratch(2)/6000;  
    Header (i,:)='linked spectrum';
    SpectraData (i,:)=fread(fp,Spectrum_Pts,'single'); % Array with individual spectra in rows
end

Wavenumbers = flip(Wavenumbers');

%Remove Data after light off
Time = Time*60; % convert to seconds

[idif indx] = min((Time - endTime).^2);
Time(indx+1:end) = [];
SpectraData(indx+1:end,:) = [];
Num_Spectra = length(Time);

%% FT-IR Assessment & plot
%% Peak 1
[ vl, lowerI ] = min( abs( Wavenumbers(1,:)-lower1 ) );
[ vu, upperI ] = min( abs( Wavenumbers(1,:)-upper1 ) );

AreaLinBase = zeros(Num_Spectra,1);
dx =length(SpectraData(1,lowerI:upperI));
for ii=1:Num_Spectra
    AreaCorrection = (SpectraData(ii,lowerI)+SpectraData(ii,upperI))*dx/2;
    AreaLinBase(ii)=sum(SpectraData(ii,lowerI:upperI))-AreaCorrection;
end
Conversion1 = 1-AreaLinBase/mean(AreaLinBase(1:80));

f1=figure;
hold off
plot(Wavenumbers,SpectraData(1,:),'LineWidth',LineThickness)
hold on
plot(Wavenumbers,SpectraData(end,:),'LineWidth',LineThickness)
xlabel(xlab2,'Fontsize',FS)
ylabel(ylab2,'Fontsize',FS)
legend('Start','End','Location','Northeast','Fontsize',FS)
colororder(Colors)
xlim([lower1 - 100 upper1 + 100])
set(gca,'fontsize',FS)
set(gca,'YTick',[])
set(gca,'XDir','reverse')
f1.Position = [1 1 figWidth figHeight];
box on
print(gcf,OutputTitleFTIR1,'-dtiff','-r1200'); 

% Add points for to show baseline
legend('AutoUpdate','off')
plot(Wavenumbers(lowerI),SpectraData(1,lowerI),'o','Color','b')
plot(Wavenumbers(lowerI),SpectraData(end,lowerI),'o','Color','b')
plot(Wavenumbers(upperI),SpectraData(1,upperI),'o','Color','b')
plot(Wavenumbers(upperI),SpectraData(end,upperI),'o','Color','b')


%% Peak 2
[ vl, lowerI ] = min( abs( Wavenumbers(1,:)-lower2 ) );
[ vu, upperI ] = min( abs( Wavenumbers(1,:)-upper2 ) );

AreaLinBase = zeros(Num_Spectra,1);
dx =length(SpectraData(1,lowerI:upperI));
for ii=1:Num_Spectra
    AreaCorrection = (SpectraData(ii,lowerI)+SpectraData(ii,upperI))*dx/2;
    AreaLinBase(ii)=sum(SpectraData(ii,lowerI:upperI))-AreaCorrection;
end
Conversion2 = 1-AreaLinBase/mean(AreaLinBase(1:80));

f2=figure;
plot(Wavenumbers,SpectraData(1,:),'LineWidth',LineThickness)
hold on
plot(Wavenumbers,SpectraData(end,:),'LineWidth',LineThickness)
xlabel(xlab2,'Fontsize',FS)
ylabel(ylab2,'Fontsize',FS)
legend('Start','End','Location','Northeast','Fontsize',FS)
colororder(Colors)
xlim([lower2 - 100 upper2 + 100])
set(gca,'fontsize',FS)
set(gca,'YTick',[])
set(gca,'XDir','reverse')
f2.Position = [1 1 figWidth figHeight];
box on
print(gcf,OutputTitleFTIR2,'-dtiff','-r1200'); 

legend('AutoUpdate','off')
plot(Wavenumbers(lowerI),SpectraData(1,lowerI),'o','Color','b')
plot(Wavenumbers(lowerI),SpectraData(end,lowerI),'o','Color','b')
plot(Wavenumbers(upperI),SpectraData(1,upperI),'o','Color','b')
plot(Wavenumbers(upperI),SpectraData(end,upperI),'o','Color','b')

%% Peak 3
[ vl, lowerI ] = min( abs( Wavenumbers(1,:)-lower3 ) );
[ vu, upperI ] = min( abs( Wavenumbers(1,:)-upper3 ) );

AreaLinBase = zeros(Num_Spectra,1);
dx =length(SpectraData(1,lowerI:upperI));
for ii=1:Num_Spectra
    AreaCorrection = (SpectraData(ii,lowerI)+SpectraData(ii,upperI))*dx/2;
    AreaLinBase(ii)=sum(SpectraData(ii,lowerI:upperI))-AreaCorrection;
end
Conversion3 = 1-AreaLinBase/mean(AreaLinBase(1:80));

%Plot
f3=figure;
plot(Wavenumbers,SpectraData(1,:),'LineWidth',LineThickness)
hold on
plot(Wavenumbers,SpectraData(end,:),'LineWidth',LineThickness)
xlabel(xlab2,'Fontsize',FS)
ylabel(ylab2,'Fontsize',FS)
legend('Start','End','Location','Northeast','Fontsize',FS)
colororder(Colors)
xlim([lower3 - 100 upper3 + 100])
set(gca,'fontsize',FS)
set(gca,'YTick',[])
set(gca,'XDir','reverse')
f3.Position = [1 1 figWidth figHeight];
box on
print(gcf,OutputTitleFTIR3,'-dtiff','-r1200'); 


legend('AutoUpdate','off')
plot(Wavenumbers(lowerI),SpectraData(1,lowerI),'o','Color','b')
plot(Wavenumbers(lowerI),SpectraData(end,lowerI),'o','Color','b')
plot(Wavenumbers(upperI),SpectraData(1,upperI),'o','Color','b')
plot(Wavenumbers(upperI),SpectraData(end,upperI),'o','Color','b')

% Kinetics Plots
Conversion = [Conversion1 Conversion2 Conversion3];
OutputTitleKinetics = {OutputTitleKinetics1 OutputTitleKinetics2 OutputTitleKinetics3};
for ii = 1:3 
    f4= figure;
    clf
plot(Time,Conversion(:,ii),'LineWidth',LineThickness)
hold on
xlabel(xlab,'Fontsize',FS)
ylabel('Acrylate Conversion','Fontsize',FS)
xlim([xlow endTime])
ylim([0 1])
set(gca,'fontsize',FS)
colororder(Colors)
% Add light
OnState = [0 1];
light = [LightOnTime LightOnTime];
plot(light,OnState,'LineWidth',3)
legend(legendLabel,'fontsize',FS,'Location',LegendPos)
f4.Position = [1 1+figHeight figWidth figHeight];
box on
print(gcf,char(OutputTitleKinetics(ii)),'-dtiff','-r1200'); 
end

%Save Data
TableOutput = [Time Conversion];
writematrix(TableOutput,OutputTitle1)
writematrix(Time,OutputTitleRawTime)
TableOutput2 = [Wavenumbers; SpectraData];
writematrix(TableOutput2,OutputTitleRaw)
