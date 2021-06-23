clc
clear all
close all
%clear all

Sample=[10 3];
File_name="D:\Université\Matrise\Article\DataOutput\symetrique\SymetricLong.xlsx";

if ~exist('data','var')
   disp('Reading data file...')
   data=xlsread(File_name);
   disp(['File has ' num2str(length(data)) ' points']);
end

%colors for phi
div_max=90;
div=0:(10/1):div_max;
n_div=length(div)-1;
color=jet(n_div);

%colors for theta
div2=[45:(10/1):135];
n_div2=length(div2)-1;
color1=jet(n_div2);

temp=data(:,5);
avg=mean(temp(temp<100)/2);
dev=std(temp(temp<100)/2);
clear phi theta
j=1;
percent=10;
reject=0;
%imax=1000;
imax=length(data);
for i=1:imax
    if i/imax>=percent/100
        disp([num2str(percent) ' %']);
        percent=percent+10;
    end

    a=data(i,4)/2; % horizontal radius (a>b)
    b=data(i,5)/2; % vertical radius
    alpha=data(i,7);
	
    if b>avg-dev && b<avg+dev && a~=b
        x0(j)=data(i,2); % x0,y0 = ellipse centre coordinates
        y0(j)=-data(i,3);
        a0(j)=a;
        b0(j)=b;
        alpha2(j)=wrapTo360(alpha-90);
        if alpha2(j)>90 && alpha2(j)<=270
            alpha2(j)=alpha2(j)-180; 
        elseif alpha>270
            alpha2(j)=alpha2(j)-360;
        end
        alpha2(j);
        beta(j)=acosd(b/a);

        r=45;
        MX=[1 0 0 ; 0 cosd(r) -sind(r); 0 sind(r) cosd(r)];
        V1=[sind(alpha2(j))*sind(beta(j)); cosd(alpha2(j))*sind(beta(j)); cosd(beta(j))];
        V2=MX*V1;

        phi(j)=acosd(V2(3));
        theta(j)=acosd(V2(2));
        z1(j,:)=color(find(div>phi(j),1)-1,:);
        z2(j,:)=color1(find(div2>theta(j),1)-1,:);
        
        % plot fibre
        plot1=0;
        if plot1==1
            figure(1)
            %plot3([0,V1(1)],[0,V1(2)],[0,V1(3)],'LineWidth',8)
            xlabel('x'),ylabel('y'),zlabel('z')
            axis equal
            hold on
            plot3([0,V2(1)],[0,V2(2)],[0,V2(3)],'r','LineWidth',5)
            view(-8.8600,-64.1000)
            locs = axis; % get current axis boundaries
            hold on; 
            plot3([locs(1) locs(2)], [0 0], [0 0]); %plot xaxis, a line between(-x,0,0) and (x,0,0);
            plot3([0 0], [locs(3) locs(4)], [0 0]); %plot y axis, the line (0,-y,0) and (0,y,0);
            plot3([0 0], [0 0], [locs(5) locs(6)]); % plot z axis
            hold off
        end
        
        j=j+1;
    else
        reject=reject+1;
    end
end

100*reject/imax

figure(1)
subplot(1,2,1)
scatter(x0,y0,10,z1,'filled')
colorbar
colormap(jet(n_div))
hcb=colorbar;
set(hcb,'YTick',div)
axis equal
axis([0 max(x0) round(min(y0)) 0])
caxis([min(div) max(div)])
title('Phi');

subplot(1,2,2)
scatter(x0,y0,10,z2,'filled')
colorbar
colormap(jet(n_div2))
hcb=colorbar;
set(hcb,'YTick',div2)
axis equal
axis([0 max(x0) round(min(y0)) 0])
caxis([min(div2) max(div2)])
title('Theta');



