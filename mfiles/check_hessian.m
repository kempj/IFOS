close all
clear all

% This m-file is for making movies of wave propagation.

%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%
% Here the basic name of the binary snapshot files must be given:
% (The default extension is *.bin. The increasing number of the snapshot
% is added to the basic name, e.g. if ten snapshots were computed
% and the basic name is snap, then the filenames for the snapshots
% are snap1.bin, snap2.bin, ..., snap10.bin.
%filediv='snap/Khang_hom.bin.p.00';

% Receiver Positions
%yrec=3.6:0.66:69.6;
%xrec=45.5.*ones(length(yrec),1);
clear all
close all

% plot distance of sources and receivers
dnrec=10;
dnsour=1;

jseafloor=440;

% write model output files
WRITEMODE=1;
READ=1;

xs = 10.0;
ys = 0.48;

%rec=load('receiver_spheres.dat');
%xrec1=rec(:,1);
%yrec1=rec(:,2);

%xrec = xrec(1:dnrec:length(xrec));
%yrec = yrec(1:dnrec:length(yrec));

%source=load('sources_plot.dat')
%xshot=source(:,1);
%yshot=source(:,3);

%xshot = xshot(1:dnsour:length(xshot));
%yshot = yshot(1:dnsour:length(yshot));

% -------------------------------------------------------------------------
% P-Wave Velocity
% -------------------------------------------------------------------------

SMOOTH=0;
IMP=1;

vpmod = [1500.0 2000.0 3000.0];
vsmod = [866.0 1155.0 1732.0]; 
rhomod = [1929.0 2073.0 2294.0]

IDX=1; IDY=1;
IDXI=10; IDYI=10;
dh=0.1;

% gridsize and grid spacing (as specified in parameter-file) 

NX1=1; NX2=500;
NY1=1; NY2=150; 

%caxis_value1 = -1e-1;
%caxis_value2 = 1e-1;

caxis_value_vp1 = -1e-15;
caxis_value_vp2 = 1e-15;

%caxis_value_vp1 = -1e-2;
%caxis_value_vp2 = 1e-2;

caxis_value_vs1 = 0;
caxis_value_vs2 = 10;

caxis_value_rho1 = -1e-1;			
caxis_value_rho2 = 1e-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grid size
nx=NX2-NX1+1 
ny=NY2-NY1+1


if(READ==1)
clear vp;
clear vs;
clear rho;

% load model
 file='jacobian_Test_p.old';
 disp([' loading file ' file]);
 fid=fopen(file,'r','ieee-le');
 grad=fread(fid,[ny,nx],'float');
 fclose(fid);

% load model
 file='jacobian_Test_hessian';
 disp([' loading file ' file]);
 fid=fopen(file,'r','ieee-le');
 hess=fread(fid,[ny,nx],'float');
 fclose(fid);

end

%vp1 = vp(1:10:ny,1:10:nx);
%clear vp;
%vp = vp1;
%clear vp1;

%vs1 = vs(1:10:ny,1:10:nx);
%clear vs;
%vs = vs1;
%clear vs1;

%rho1 = rho(1:10:ny,1:10:nx);
%clear rho;
%rho = rho1;
%clear rho1;

NX1=1; NX2=500;
NY1=1; NY2=150; 

dh=20.0;

% grid size
nx=NX2-NX1+1 
ny=NY2-NY1+1

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;

epsilon = 5;
maxhess = max(max(hess));
  for i=1:nx
      for j=1:ny
          
          depth = j.*dh;
          ihess(j,i) = depth/(hess(j,i)+epsilon);
          
      end
  end


% apply damping of H^-1 in the PML layers
DAMPING=99;
amp=1.0-DAMPING./100.0;   % amplitude at the edge of the numerical grid 
ifw=10;  % frame width in gridpoints 
a=sqrt(-log(amp)./(ifw.^2));
        
for i=1:ifw
    coeff(i)=exp(-(a.^2.*(ifw-i).^2));
    %coeff(i) = 0.0;
end

% initialize array of coefficients with one 
for j=1:ny
    for i=1:nx 
      absorb_coeff(j,i)=1.0;
    end
end

% coefficients for left and right boundary
yb=1; 
ye=ny;
  for i=1:ifw
      ye=ny-i+1;
        for j=yb:ye
          absorb_coeff(j,i)=coeff(i);
        end           
  end

yb=1; 
ye=ny;
  for i=1:ifw
      ii=nx-i+1;
      ye=ny-i+1;
        for j=yb:ye
          absorb_coeff(j,ii)=coeff(i);
        end           
  end
  
% coefficients for bottom boundary
xb=1; 
xe=nx;
  for j=1:ifw
      xb=j;
      jj=ny-j+1;
      xe=nx-j+1;
        for i=xb:xe
          absorb_coeff(jj,i)=coeff(j);
        end           
  end

 % apply damping in PML boundaries 
  ihess = ihess .* absorb_coeff;

    subplot 311
	imagesc(x,y,grad);
	
	hold on;
    plot(xs,ys,'k+');	
       %caxis([0 1]);
    %caxis([caxis_value_vp1 caxis_value_vp2]);
        %colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       %set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
		 axis ij
		 axis equal
		 axis tight
       
       ylabel('y [m]');
       
       iter_text=['Gradient'];
       
       title(iter_text);


subplot 312
imagesc(x,y,ihess);
	
	hold on;
plot(xs,ys,'k+');	
       %caxis([caxis_value_vs1 caxis_value_vs2]);
    
        %colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       %set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
		 axis ij
                 axis equal
		 axis tight
		 

       ylabel('y [m]');
       
       iter_text=['P (Hessian + 5I)^{-1}'];
       
       title(iter_text);

%if(WRITEMODE==1)       
%file1=['sphere_complex.vs'];
%fid1=fopen(file1,'w','ieee-le');
%fwrite(fid1,vs,'float')
%fclose(fid1);
%end

subplot 313
imagesc(x,y,ihess.*grad);
	
	hold on;
plot(xs,ys,'k+');	
       % caxis([caxis_value_rho1 caxis_value_rho2]);
    
        %colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       %set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
	
                 axis equal
		 axis tight
       xlabel('x [m]');
       ylabel('y [m]');
       
       iter_text=['P (H + 5 I)^{-1} * Gradient'];
       
       title(iter_text);

if(WRITEMODE==1)       
file1=['taper.bin'];
fid1=fopen(file1,'w','ieee-le');
fwrite(fid1,ihess,'float')
fclose(fid1);
end       

figure;
imagesc(x,y,ihess);
	
	 hold on;
     plot(xs,ys,'k+');	
     % caxis([caxis_value_rho1 caxis_value_rho2]);
    
        colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       %set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
	
                 axis equal
		 axis tight
       xlabel('x [m]');
       ylabel('y [m]');
       
       iter_text=['H^{-1}'];
       
       title(iter_text);

figure;
imagesc(x,y,absorb_coeff);
	
	 hold on;
     plot(xs,ys,'k+');	
     % caxis([caxis_value_rho1 caxis_value_rho2]);
    
        colorbar;
	%set(gca,'YDir','normal');
        %colorbar;
        %axis equal;
       %set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
	
                 axis equal
		 axis tight
       xlabel('x [m]');
       ylabel('y [m]');
       
       iter_text=['damping coeff'];
       
       title(iter_text);

%if(WRITEMODE==1)       
%file1=['sphere_complex.rho'];
%fid1=fopen(file1,'w','ieee-le');
%fwrite(fid1,rho,'float')
%fclose(fid1);
%end
