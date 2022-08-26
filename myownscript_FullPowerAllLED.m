theta = 70;
% semi-angle at half power
m=-log10(2)/log10(cosd(theta));
%Lambertian order of emission
P_LED=20;
%transmitted optical power by individual LED
nLED=60;
% number of LED array nLED*nLED 
P_total=nLED*nLED*P_LED;
%Total transmitted power 
Adet=1e-4;
%detector physical area of a PD
Ts=1;
%gain of an optical filter; ignore if no filter is used
index=1.5;
%refractive index of a lens at a PD; ignore if no lens is used
FOV=60;
%FOV of a receiver
G_Con=(index^2)/(sind(FOV).^2);
%gain of an optical concentrator; ignore if no lens is used
%%
lx=5; ly=5; lz=3;
% room dimension in meter
h=1.65;
%the distance between source and receiver plane
%[XT,YT]=meshgrid([-lx/4 lx/4],[-ly/4 ly/4]);
[LED_OrgX,LED_OrgY] = meshgrid(-0.3:0.01:0.29);
[LED_A_X,LED_A_Y] = deal(LED_OrgX-lx/4,LED_OrgY-ly/4);
[LED_B_X,LED_B_Y] = deal(LED_OrgX+lx/4,LED_OrgY-ly/4);
[LED_C_X,LED_C_Y] = deal(LED_OrgX-lx/4,LED_OrgY+ly/4);
[LED_D_X,LED_D_Y] = deal(LED_OrgX+lx/4,LED_OrgY+ly/4);
% position of LED; it is assumed all LEDs are located at same point for
% faster simulation
% for one LED simulation located at the central of the room, use XT=0 and YT=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=lx*5; Ny=ly*5;
% number of grid in the receiver plane
x=linspace(-lx/2,lx/2,Nx); y=linspace(-ly/2,ly/2,Ny); [XR,YR]=meshgrid(x,y);
% grid of receivers
P_rec_A = zeros(size(XR));
% matrix for accumulating receiver power.
for LED_Row = 1:nLED
    for LED_Colum = 1:nLED
        D1=sqrt((XR-LED_A_X(LED_Row,LED_Colum)).^2+(YR-LED_A_Y(LED_Row,LED_Colum)).^2+h^2);
        % distance vector from source 1
        cosphi_A1 = h./D1;
        % angle vector
        receiver_angle = acosd(cosphi_A1);
        % alternative methods to calculate angle, more accurate if theangle are
        % negatives
        H_A1=(m+1)*Adet.*cosphi_A1.^(m+1)./(2*pi.*D1.^2);
        % channel DC gain for source 1
        P_rec_A1=P_LED.*H_A1.*Ts.*G_Con;
        % received power from source 1;
        P_rec_A1(find(abs(receiver_angle)>FOV))=0;
        % if the anlge of arrival is greater than FOV, no current is generated at
        % the photodiode.
        P_rec_A = P_rec_A + P_rec_A1;
        % acculumate received power from every LED in Zone A
    end
end

P_rec_B=fliplr(P_rec_A);
% received power from LED Zone B, due to symmetry no need separate % calculations
P_rec_C=flipud(P_rec_A);
P_rec_D=fliplr(P_rec_C);
P_rec_total=P_rec_A+P_rec_B+P_rec_C+P_rec_D;



rho=0.8;
%reflection coefficient


lx=5; ly=5; lz=1.65;
% room dimension in meter
%[XT,YT,ZT]=meshgrid([-lx/4 lx/4],[-ly/4 ly/4],lz/2); % position of Transmitter (LED);

[LED_OrgX,LED_OrgY,LED_OrgZ] = meshgrid(-0.3:0.01:0.29,-0.3:0.01:0.29,lz/2);
[LED_A_X,LED_A_Y,LED_A_Z] = deal(LED_OrgX-lx/4,LED_OrgY-ly/4,LED_OrgZ);
[LED_B_X,LED_B_Y,LED_B_Z] = deal(LED_OrgX+lx/4,LED_OrgY-ly/4,LED_OrgZ);
[LED_C_X,LED_C_Y,LED_C_Z] = deal(LED_OrgX-lx/4,LED_OrgY+ly/4,LED_OrgZ);
[LED_D_X,LED_D_Y,LED_D_Z] = deal(LED_OrgX+lx/4,LED_OrgY+ly/4,LED_OrgZ);

Nx=lx*5; Ny=ly*5; Nz=round(lz*5);
% number of grid in each surface
dA=lz*ly/(Ny*Nz);
% calculation grid area
x=linspace(-lx/2,lx/2,Nx); y=linspace(-ly/2,ly/2,Ny); z=linspace(-lz/2,lz/2,Nz); [XR,YR,ZR]=meshgrid(x,y, -lz/2);
%%
%first transmitter calculation
%TP1=[XT(1,1,1) YT(1,1,1) ZT(1,1,1)];
TP1=[LED_A_X(1,1,1) LED_A_Y(1,1,1) LED_A_Z(1,1,1)];
h11=zeros(Nx,Ny);
h22=zeros(Nx,Ny);
h33=zeros(Nx,Ny);
h44=zeros(Nx,Ny);
h1=zeros(Nx,Ny);
h2=zeros(Nx,Ny);
h3=zeros(Nx,Ny);
h4=zeros(Nx,Ny);

for mm = 1:60
    disp(mm);
    for nn = 1:60
        TP1=[LED_A_X(mm,nn,1) LED_A_Y(mm,nn,1) LED_A_Z(1,1,1)];
        
        % transmitter position
        TPV=[0 0 -1];
        
        % transmitter position vector
        RPV=[0 0 1];
        % receiver position vector
        %%
        %%%%%%%%%%%%%%%calculation for wall 1%%%%%%%%%%%%%%%%%%
        WPV1=[1 0 0];
        % position vector for wall 1
        for ii=1:Nx
          for jj=1:Ny
             RP=[x(ii) y(jj) -lz/2];
             % receiver position vector
             h1(ii,jj)=0;
             % reflection from North face
             for kk=1:Ny
                for ll=1:Nz
                    WP1=[-lx/2 y(kk) z(ll)]; 
                    D1=sqrt(dot(TP1-WP1,TP1-WP1)); 
                    cos_phi= abs(WP1(3)-TP1(3))/D1; 
                    cos_alpha=abs(TP1(1)- WP1(1))/D1;
                    D2=sqrt(dot(WP1-RP,WP1-RP)); 
                    cos_beta=abs(WP1(1)- RP(1))/D2; 
                    cos_psi=abs(WP1(3)- RP(3))/D2;
                    if abs(acosd(cos_psi))<=FOV 
                        h1(ii,jj)=h1(ii,jj)+(m+1)*Adet*rho*dA*cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                    end 
                end
             end 
          end
        end
        h11 = h11 + h1;
        %%
        WPV2=[0 1 0];
        % position vector for wall 2
        %% calculation of the channel gain is similar to wall1 
        for ii=1:Nx
          for jj=1:Ny
             RP=[x(ii) y(jj) -lz/2];
             % receiver position vector
             h2(ii,jj)=0;
             % reflection from wall 2
             for kk=1:Ny
                for ll=1:Nz
                    WP1=[lx/2 y(kk) z(ll)]; 
                    D1=sqrt(dot(TP1-WP1,TP1-WP1)); 
                    cos_phi= abs(WP1(3)-TP1(3))/D1; 
                    cos_alpha=abs(TP1(1)- WP1(1))/D1;
                    D2=sqrt(dot(WP1-RP,WP1-RP)); 
                    cos_beta=abs(WP1(1)- RP(1))/D2; 
                    cos_psi=abs(WP1(3)- RP(3))/D2;
                    if abs(acosd(cos_psi))<=FOV 
                        h2(ii,jj)=h2(ii,jj)+(m+1)*Adet*rho*dA*cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                    end 
                end
             end 
          end
        end
        h22 = h22 + h2;
        
        for ii=1:Nx
          for jj=1:Ny
             RP=[x(ii) y(jj) -lz/2];
             % receiver position vector
             h3(ii,jj)=0;
             % reflection from North face
             for kk=1:Ny
                for ll=1:Nz
                    WP1=[x(kk) -ly/2 z(ll)]; 
                    D1=sqrt(dot(TP1-WP1,TP1-WP1)); 
                    cos_phi= abs(WP1(3)-TP1(3))/D1; 
                    cos_alpha=abs(TP1(1)- WP1(1))/D1;
                    D2=sqrt(dot(WP1-RP,WP1-RP)); 
                    cos_beta=abs(WP1(1)- RP(1))/D2; 
                    cos_psi=abs(WP1(3)- RP(3))/D2;
                    if abs(acosd(cos_psi))<=FOV 
                        h3(ii,jj)=h3(ii,jj)+(m+1)*Adet*rho*dA*cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                    end 
                end
             end 
          end
        end
        h33 = h33 + h3;
        
        for ii=1:Nx
          for jj=1:Ny
             RP=[x(ii) y(jj) -lz/2];
             % receiver position vector
             h4(ii,jj)=0;
             % reflection from North face
             for kk=1:Ny
                for ll=1:Nz
                    WP1=[x(kk) ly/2 z(ll)]; 
                    D1=sqrt(dot(TP1-WP1,TP1-WP1)); 
                    cos_phi= abs(WP1(3)-TP1(3))/D1; 
                    cos_alpha=abs(TP1(1)- WP1(1))/D1;
                    D2=sqrt(dot(WP1-RP,WP1-RP)); 
                    cos_beta=abs(WP1(1)- RP(1))/D2; 
                    cos_psi=abs(WP1(3)- RP(3))/D2;
                    if abs(acosd(cos_psi))<=FOV 
                        h4(ii,jj)=h4(ii,jj)+(m+1)*Adet*rho*dA*cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                    end 
                end
             end 
          end
        end
        h44 = h44 + h4;

    end
end

P_rec_A1=(h11+h22+h33+h44)*P_LED.*Ts.*G_Con;
% h2, h3 and h4 are channel gain for walls 2,3 and 4, respectively.
P_rec_A2=fliplr(P_rec_A1);
% % received power from source 2, due to symmetry no need separate % calculations
P_rec_A3=flipud(P_rec_A1);
P_rec_A4=fliplr(P_rec_A3);
P_rec_total_1ref=P_rec_A1+P_rec_A2+P_rec_A3+P_rec_A4;



P_rec_total = P_rec_total + P_rec_total_1ref;


P_rec_dBm=10*log10(P_rec_total);
%%%%%%%%%%%%%%%%%%%%
figure;
surfc(x,y,P_rec_dBm);
colorbar;
xlabel('Length of room [m]')
ylabel('Width of room [m]')
zlabel('Received Power in (dBm)')