% AUTHOR: Hugues Berry
%         hugues.berry@inria.fr
%         http://www.inrialpes.fr/Berry/
%
% REFERENCE: Blum Moyse & Berry. Modelling the modulation of cortical Up-Down state switching by astrocytes
%
% LICENSE: CC0 1.0 Universal

clear all; close all;
tmax=6e3;
dt=0.2;
tt=(0:(tmax/dt))*dt;

newstart=1;
thresh=50;%1.25;%thresh to separate UP from Down states
U = []; D = []; PER_UP = []; PER_UP2 = []; PER_DO2 = [];

for numm=0:1:199
   
    %folder = './data/Jai0_outrate/';
    folder = './data/Jai_sup_0_outrate/';
   
    Y=csvread(char(strcat(folder,'outRateBuff',num2str(numm),'.csv')));
    
    [per_UP,u,d] = percentUP(Y(newstart:end),thresh,numm);
    U = [U,reshape(u,[1,length(u)])];
    D = [D,reshape(d,[1,length(d)])];
    
end
'# of U'
numel(U)
'# of D'
numel(D)
'mean perup'
mean(U)
'std perup'
std(U)
'mean perdo'
mean(D)
'std perdo'
std(D)


h = figure(1);
h.Position= [0 0 1500 400];
histogram(U/100,linspace(0,max(U/100),100),'FaceColor','#EDB120');
hold on;
histogram(D/100,linspace(0,max(U/100),100),'FaceColor','#7E2F8E');
legend();
ax = gca;
ax.FontSize = 20;
xlim([0,8]);
xlabel('duration (s)','Fontsize',20)
ylabel('no. event', 'Fontsize',20)
hold off;

function [per_UP,U,D]=percentUP(data,thresh,numm,sim,simmax)

cross=0;
crosses=[];
smootheddata=smoothdata(data,'movmedian',10);
currdata=smootheddata;

while ~isempty(cross)
    if currdata(1)>thresh
        cross=-1+find((currdata<thresh)>0,1);
    else
        cross=-1+find((currdata>thresh)>0,1);
    end
    
    crosses=[crosses, cross];
    currdata=currdata(cross+1:end);
end
swaps=[1 cumsum(crosses)];
bound_c1=[];
bound_c2=[];
val_c1=[];
val_c2=[];
for i=1:2:(numel(swaps)-1)
   bound_c1=[bound_c1 [swaps(i) swaps(i+1)-1]];
   if i<=numel(swaps)-2
       bound_c2=[bound_c2 [swaps(i+1) swaps(i+2)-1]];
   else
       bound_c2=[bound_c2 [swaps(i+1)]];
   end
end

if rem(numel(swaps),2)==0
    bound_c2=[bound_c2 numel(data)];
else
    bound_c1=[bound_c1 swaps(end) numel(data)];
end


for i=1:2:numel(bound_c1)-1
   val_c1=[val_c1 data(bound_c1(i):bound_c1(i+1))];
end
for i=1:2:numel(bound_c2)-1
   val_c2=[val_c2 data(bound_c2(i):bound_c2(i+1))];
end

av_c1=mean(val_c1);
av_c2=mean(val_c2);

hold on
time_c1=0;
time_c2=0;
Ltime_c1 = zeros(1,length(bound_c1));
Ltime_c2 = zeros(1,length(bound_c2));
for i=1:2:numel(bound_c1)
    time_c1=time_c1+bound_c1(i+1)-bound_c1(i)+1;
    Ltime_c1(i) = bound_c1(i+1)-bound_c1(i)+1;
end
for i=1:2:numel(bound_c2)
    time_c2=time_c2+bound_c2(i+1)-bound_c2(i)+1;
    Ltime_c2(i) = bound_c2(i+1)-bound_c2(i)+1;
end

if av_c1>av_c2
    per_UP=time_c1/numel(data)*100;  
    U = nonzeros(Ltime_c1(2:end-2));
    D = nonzeros(Ltime_c2(2:end-2));

else
    per_UP=time_c2/numel(data)*100;
    U = nonzeros(Ltime_c2(2:end-2));
    D = nonzeros(Ltime_c1(2:end-2));
end
if per_UP==0 && ((av_c1>3*thresh)||(av_c2>3*thresh))%only the UP phase
   per_UP=100; 
end
if (nargin>2)
    disp(['time in up=',num2str(per_UP),'%']);
end

% figure();
% plot(smootheddata);
% hold on;
% %plot(bound_c1,zeros(1,length(bound_c1)),'r*');
% plot(bound_c2,zeros(1,length(bound_c2)),'b*');
% hold off;
% title(num2str(numm))
end

