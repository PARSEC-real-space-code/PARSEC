function [dd]=bands(figid,ofs)

% this script is intentended for plotting the band structure
% results of PARSEC calculations that were written to the file bands.dat 
% the script was written by Noa Marom and modified by Amir Natan. 

%input parameters:

%figid - figure number in matlab. should be a positive integer. 
%ofs - offset in eV if desired, can be usually 0.0 .
%
% a simple call should look like dd=bands(1,0);

%output parameters:

%dd - an array containing the band data.

%filename= name of parsec bands output file including extention
%nbands= number of states from the file parsec.in
%units - Rydberg, Hartree or eV


filename = input('Enter file name: ','s');
nbands = input('Enter number of bands: ');
units = input('Select units- enter Ry, Ha , or eV: ', 's'); 
plot_color = input('Select color scheme: 0 for b/w 1 for rainbow: ');

figure(figid);
last=6+nbands;
fid=fopen(filename,'r');

i=1;
while(1)
[d1,count]=fscanf(fid,'%f',last);
if(count==last)
   dd(i,:)=d1;
   i=i+1;
else
    break
end;
end;

switch (units)
    case 'Ha'
        dd(:,7:last)=dd(:,7:last)/2;
    case 'eV'
        dd(:,7:last)=dd(:,7:last)*13.6055+ofs;
end

switch(plot_color)
    case 1
        % color order is: purple navy blue cyan green yellow orange
        %red magenta maroon
        colors=[0.5  0.0  0.5
                0.0  0.0  0.5
                0.0  0.0  1.0
                0.0  1.0  1.0
                0.0  1.0  0.0
                1.0  1.0  0.0
                1.0  0.75  0.0
                1.0  0.0  0.0
                1.0  0.0  1.0
                0.5  0.0  0.0];
        
        set(gcf,'DefaultAxesColorOrder',colors)
        plot(dd(:,7:last),'LineWidth',2);
        hold on;
    case 0
        figure(figid);
        plot(dd(:,7:last),'k','LineWidth',2);
        hold on;
end

dum = ylim;
dum1=[dum(1):0.01:dum(2)];
dum2=size(dd);

for ii=1:(dum2(1)-1);
    if (dd(ii+1,2) ~= dd(ii,2))
        dum3 = ones(size(dum1))*(ii+1);
        figure(figid); plot(dum3,dum1,'k');
    end
end

hold off;


