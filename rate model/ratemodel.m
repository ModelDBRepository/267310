% AUTHOR: Hugues Berry
%         hugues.berry@inria.fr
%         http://www.inrialpes.fr/Berry/
%
% REFERENCE: Blum Moyse & Berry. Modelling the modulation of cortical Up-Down state switching by astrocytes
%
% LICENSE: CC0 1.0 Universal

%%% Fixed Parameters %%
% time in msec
taue=10;
taui=2;
taua=20;
tauadap=500;
thetai=25;
thetaa=-3.5;
Jee=5;
Jei=-1;
Jii=-0.5;
Jie=10;
Jaa=0.1;
Jea=1;
Jia=0.5;
ge=1;
gi=4;
ga=1;
%%%%%%%
thetae=20;
beta=2;
Jae=0.5;
Jai=0.5;

%%main composite parameters %%
Jpii=Jii-1/gi;
Jpaa=Jaa-1/ga;
Jpee=Jee-1/ge;
ha=Jpee*Jpii-Jei*Jie;
he=Jpii*Jpaa-Jia*Jai;
hi=Jpaa*Jpee-Jae*Jea;
fia=Jie*Jea-Jia*Jpee;
fei=Jea*Jai-Jei*Jpaa;
fea=Jei*Jia-Jea*Jpii;
fie=Jia*Jae-Jie*Jpaa;
fae=Jai*Jie-Jae*Jpii;
fai=Jae*Jei-Jai*Jpee;

%weight matrix
Jmat=[Jee Jei Jea;Jie Jii Jia;Jae Jai Jaa];

%paramter vectors
tauvec=[taue;taui;taua];
gvec=[ge;gi;ga];

%%param for the runge kutta ODE integration
tmax=6e3;
dt=0.2;
tt=(0:(tmax/dt))*dt;

%number of simulations per point (to average % time in Up state)
nsims=10;
%remove the first eight of the simulation to get stationnary dynamics
newstart=round((ceil(tmax/dt)+1)/8);

%thresh to separate UP from Down states based on rI
thresh=1.25;

%parameter ranges to be explored in terms of thetaE and beta values
allthetaes=-10:0.5:20;
allbetas=0:0.5:10;

%the matrix that contains de % of time in Up state for each (thetaE,beta)
%couples explored
allperUP=zeros(numel(allthetaes),numel(allbetas));

for i=1:numel(allthetaes)
    thetae=allthetaes(i)
    for j=1:numel(allbetas)
        beta=allbetas(j)
        thetavec=[thetae;thetai;thetaa];
        paradap=[tauadap,beta];
        per_UP=0;
        for k=1:nsims
            % simulate nsims time the ODE system with the current thetaE
            % and beta and random initial conditions take as random
            % integers between 0 and 4
            Y=rateRK(randi(5,1,4)-1,tmax,dt,Jmat,thetavec,gvec,tauvec,paradap);
            %comput the percentage of time spent in the UP state for the
            %r_I time-series of the current simulation
            per_UP=per_UP+percentUP(Y(2,newstart:end),thresh);
        end
        %average over nsims and save in allperUP
        per_UP=per_UP/nsims;
        allperUP(i,j)=per_UP;
        disp(['% in UP=',num2str(per_UP)]);
    end
end

%plot the color map
imagesc(allthetaes,allbetas,allperUP');ylabel('$\beta$','Interpreter','latex');xlabel('$\theta_E$','Interpreter','latex');title('% time in UP');
colormap(pink);colorbar;caxis([0 100]);
set(gca,'YDir','normal');

%compute the theoretical boundaries between the different regimes
eq11=-ga*Jea*thetaa/(1-ga*Jaa);
eq20=-(hi*thetai+fia*thetaa)/fie;
eq19=(fie*allthetaes+hi*thetai+fia*thetaa)/(Jpaa*thetai-Jia*thetaa);
eq21=-allthetaes*(Jae*fea+Jai*fia+Jpaa*ha)./(fei*thetai+fea*thetaa+allthetaes*(he-Jpaa*Jpii-Jai*Jia));

%plot the boundaries
hold on;
plot([eq11 eq11],[min(allbetas) max(allbetas)], '--','Color','#00FFFF','LineWidth',1);
plot([eq20 eq20],[min(allbetas) max(allbetas)], '--','Color','#77AC30','LineWidth',1);
plot(allthetaes,eq19, '-.','Color','#FF00FF','LineWidth',1);
plot(allthetaes,eq21, 'k--','LineWidth',1);legend('eq11','eq20','eq19','eq21');
hold off;


function per_UP=percentUP(data,thresh,sim,simmax)
%the functin that computes the time spent in the UP state of the 
%time-series given as input
cross=0;
crosses=[];

%1. smooth data using a sliding window
%nb: the fonction smoothdata is available only from version R2017a
smootheddata=smoothdata(data,'movmedian',100);
currdata=smootheddata;

while ~isempty(cross)
    if currdata(1)>thresh
        %if currdata starts in UP state
        %find the end of the UP state as the first crossing of the
        %threshold from above
        cross=-1+find((currdata<thresh)>0,1);
    else
        %if currdata starts in DOWN state
        %find the end of the DOWN state as the first crossing of the
        %threshold from below
        cross=-1+find((currdata>thresh)>0,1);
    end
    %add the crossing/transition time
    crosses=[crosses, cross];
    %remove the current phase and start with the next one
    currdata=currdata(cross+1:end);
end

% restore the initial crossing/transition times
swaps=[1 cumsum(crosses)];
bound_c1=[];
bound_c2=[];
val_c1=[];
val_c2=[];

for i=1:2:(numel(swaps)-1)
   %left and right boundaries of the UP or DOWN phases
   bound_c1=[bound_c1 [swaps(i) swaps(i+1)-1]];
   %left and right boundaries of the DOWN or UP phases
   if i<=numel(swaps)-2
       bound_c2=[bound_c2 [swaps(i+1) swaps(i+2)-1]];
   else
       bound_c2=[bound_c2 [swaps(i+1)]];
   end
end

%correct for the parity in the number of phases
if rem(numel(swaps),2)==0
    bound_c2=[bound_c2 numel(data)];
else
    bound_c1=[bound_c1 swaps(end) numel(data)];
end


%do not take into account the very last phase (since it is not complete)
if ~isempty(bound_c2) && ~isempty(bound_c1)
    if max(bound_c2)>max(bound_c1)
        bound_c2=bound_c2(1:(end-2));
    else
        bound_c1=bound_c1(1:(end-2));
    end
end

%durations of the corresponding phases
for i=1:2:numel(bound_c1)-1
   val_c1=[val_c1 data(bound_c1(i):bound_c1(i+1))];
end
for i=1:2:numel(bound_c2)-1
   val_c2=[val_c2 data(bound_c2(i):bound_c2(i+1))];
end

%averages
av_c1=mean(val_c1);
av_c2=mean(val_c2);

hold on
time_c1=0;
time_c2=0;

for i=1:2:numel(bound_c1)
    %if plots are needed
    if (nargin>2)
        subplot(simmax,1,sim),plot([bound_c1(i) bound_c1(i+1)],[av_c1 av_c1],'linewidth',3);
    end
    time_c1=time_c1+bound_c1(i+1)-bound_c1(i)+1;
end
for i=1:2:numel(bound_c2)
    if (nargin>2)
        subplot(simmax,1,sim),plot([bound_c2(i) bound_c2(i+1)],[av_c2 av_c2],'linewidth',3);
    end
    time_c2=time_c2+bound_c2(i+1)-bound_c2(i)+1;
end
ylim([-1,inf]);
hold off

%if c1=UP and c2=DOWN
if av_c1>av_c2
    per_UP=time_c1/(time_c1+time_c2)*100;   
else %if c2=UP and c1=DOWN
    per_UP=time_c2/(time_c1+time_c2)*100;
end

% account for the fact that one may have only a unique UP phase
if per_UP==0 && ((av_c1>2*thresh)||(av_c2>2*thresh))%only the UP phase
   per_UP=100; 
end

if (nargin>2)
    disp(['time in up=',num2str(per_UP),'%']);
end
end


function y=rateRK(rinit,tmax,dt,Jmat,thetavec,gvec,tauvec,paradap)
%Runge-Kutta numerical resultion of the rate model
%classical RK4 scheme for Stratonovitch SDEs 
%see eg https://en.wikipedia.org/wiki/Runge–Kutta_method_(SDE)


hstep=dt/2;
nsteps=tmax/dt;
y=zeros(4,nsteps+1);
y(:,1)=rinit;

%noise
sigman=3.5*sqrt(2);
sigman2=sigman^2;
taun=1;
noise=zeros(3,1);
prefn1=exp(-dt/taun);
prefn2=sqrt(sigman2*taun/2*(1-exp(-2*dt/taun)));
checknoise=zeros(3,nsteps);

for i=1:nsteps
    t=i*dt;
    k_1 = ratemodelODE(t,y(:,i),Jmat,thetavec,gvec,tauvec,paradap,noise);
    k_2 = ratemodelODE(t+hstep,y(:,i)+hstep*k_1,Jmat,thetavec,gvec,tauvec,paradap,noise);
    k_3 = ratemodelODE(t+hstep,y(:,i)+hstep*k_2,Jmat,thetavec,gvec,tauvec,paradap,noise);
    k_4 =ratemodelODE(t+dt,y(:,i)+dt*k_3,Jmat,thetavec,gvec,tauvec,paradap,noise);
    y(:,i+1) = y(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    noise=prefn1*noise+prefn2*randn(3,1);    
    checknoise(:,i)=noise;
end

end

function dy = ratemodelODE(t,y,Jmat,thetavec,gvec,tauvec,paradap,noise)
%the ODEs of the rate model

dy=zeros(4,1);
V=Jmat*y(1:3,1)-thetavec+noise;
V(1,1)=V(1,1)-y(4,1);

%rectification of the input
V(V<0)=0;
%ODEs for the rate functions
dy(1:3,1)=((gvec.*V)-y(1:3,1))./tauvec;
%ODE for the adaptation
dy(4,1)=(-y(4,1)+paradap(2)*y(1))/paradap(1);

end
