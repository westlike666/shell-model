%This part solves the rate equations and saves results in workspace
clear; clc;
N=50;%number of shells
t_max=200;
steps=4000;
frame_rate=25; %video framerate
global sigma_z sigma_y sigma_x d_p

sigma_z=0.42;%mm
sigma_y=sigma_z;
sigma_x=0.75;%mm
n=65; %PQN
d_p=0.3; %peak density in um-3

sigma_env=2;%consider amount of sigma environments

pos=linspace(0,sigma_env*sigma_z-0.5*sigma_env*sigma_z/(N-0.5),N);
pos_border=linspace(0.5*sigma_env*sigma_z/(N-0.5),sigma_env*sigma_z,N);
d=arrayfun(@(z) d_p*exp(-(z^2)/(2*sigma_z^2)),pos);
border=arrayfun(@(z) d_p*exp(-(z^2)/(2*sigma_z^2)),pos_border);
vol=4/3*pi*pos_border.^3*sigma_x/sigma_z*10^9;
vol(2:end)= diff(vol,[],2);
e_total=sum(d.*vol);

mkdir(date());
filename=[date(),'\','shell_sim_n=',num2str(n),'_d0=', num2str(d_p), '_sigma_', num2str(sigma_z), 'mm_', num2str(sigma_x),'mm','_tfinal',num2str(t_max),'ns_',num2str(sigma_env),'sigmaenv'];
%solve all rate equations
t1=clock;
[time,nden,eden,deac,Te,y0]=shell_rate_eqn_sim(d, vol,n, t_max/steps, t_max) ;
save(strcat([filename, '.mat']))
t2=clock;
computation_time_in_min=sum((t2-t1).*[0,0,24*60,60, 1, 1/60])

%%
%Display data
n_x=500;%resolution in x
n_y=500;%resolution in y

x=linspace(-sigma_x*sigma_env,sigma_x*sigma_env,n_x);
y=linspace(-sigma_y*sigma_env,sigma_y*sigma_env,n_y);

M=zeros(n_y,n_x);
MR=zeros(n_y,n_x);


v=VideoWriter(strcat([filename, '.avi']));
v.FrameRate=frame_rate;
open(v)

figure('position',[0,0,1150,800])
first=true;
len=length(border);
for t=0:length(time)
    if(t==0)%find max values
        t=length(time);
    end
    for i=1:n_y
        for j=1:n_x
            for k=1:len
                den=d_p*exp(-(x(j)^2)/(2*sigma_x^2)-(y(i)^2)/(2*sigma_y^2)-(0^2)/(2*sigma_z^2));
                if den>border(k) 
                    M(i,j)=eden(t,k);
                    MR(i,j)=eden(t,k)/d(k);
                    break
                elseif k==length(border)
                    break
                end
            end
        end
    end

 

    if(first==false)
        subplot(3,2,1)
        imagesc(x,y,M',lim1)
        colorbar()
        xlabel('x (mm)','Interpreter','Latex')
        ylabel('y (mm)','Interpreter','Latex')
        title(strcat(['density in $\mu m^{-3}$ for $t=',num2str(round(time(t))),'$ ns']),'Interpreter','Latex','fontsize',12);
        axis equal
        colormap jet

        subplot(3,2,2)
        imagesc(x,y,MR',lim2)
        colorbar()
        xlabel('x (mm)','Interpreter','Latex')
        ylabel('y (mm)','Interpreter','Latex')
        title(strcat(['rel. density for $t=',num2str(round(time(t))),'$ ns']),'Interpreter','Latex','fontsize',12);
        axis equal
        colormap jet

        subplot(3,2,3)
        den_shell=eden(t,:);
        bar(den_shell,1)
        xlabel('shell number','Interpreter','Latex')
        ylabel('density in $\mu m^{-3}$','Interpreter','Latex')
        ylim(lim3);

        subplot(3,2,4)
        bar(den_shell./d,1)
        xlabel('shell number','Interpreter','Latex')
        ylabel('rel. density','Interpreter','Latex')
        ylim(lim4);

        subplot(3,2,5)
        bar(den_shell.*vol,1)
        xlabel('shell number','Interpreter','Latex')
        ylabel('number of electrons $n_e$','Interpreter','Latex')
        ylim(lim5);

        subplot(3,2,6)
        bar(den_shell.*vol./e_total,1)
        xlabel('shell number','Interpreter','Latex')
        ylabel('rel. number of electrons','Interpreter','Latex')
        ylim(lim6);
    end
    if(first)
        lim1=[0 max(max(M))];
        lim2=[0 max(max(MR))];
        lim3=[0 max(eden(t,:))*1.05];
        lim4=[0 max(eden(t,:)./d)*1.05];
        lim5=[0 max(eden(t,:).*vol)*1.05];
        lim6=[0 max(eden(t,:).*vol./e_total)*1.05];
        first=false;
        t=0;
    end
     drawnow
    frame = getframe(gcf);
    writeVideo(v, frame);
end


close(v)

%%
figure('position',[0,0,650,400])
subplot(1,2,1)
plot(time,Te);
title('Temperature evolution');
xlabel('time $t$ in ns','Interpreter','Latex')
ylabel('electron temperature $T_e$','Interpreter','Latex')

subplot(1,2,2)
for l=1:N
    plot(time, eden(:,l)/d(l));
    hold on;
end
title('Avalanche for all shells');
xlabel('time $t$ in ns','Interpreter','Latex')
ylabel('relative electron density $\rho_e/\rho_0$','Interpreter','Latex')

saveas(gcf,strcat([filename,'.pdf']));
saveas(gcf,strcat([filename,'.fig']));
