function bands2(fignum)

%this function plots spin-polarized band structure from 
%the file bands.dat that was produced by PARSEC.
%the parameter fignum is the figure id of matlab - should be 1 or bigger
%integer. 
%
%the function was written by Noa Marom and modified by Amir Natan. 
%

%bands for a spin polarized system
%up is plotted in blue, down in red
%filename= name of parsec bands output file including extention
%nbands= number of eigenvalues given for each k-point in the file

filename = input('Enter file name: ','s');
nbands = input('Enter number of bands: ');
units = input('Select units- enter Ry, Ha , or eV: ', 's'); 
plot_type = input('Select plot type: 0 for superimposed 1 for separate: ');

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
        dd(:,7:last)=dd(:,7:last)*13.6055;
end

dum2=size(dd);
dum4= dum2(1)/2;
last2=dum2(1);
switch(plot_type)
    case 0
        figure(fignum); 
        plot(dd(dum4+1:last2,7:last),'r','linewidth',2);
        hold on;
        plot(dd(1:dum4,7:last),'b','linewidth',2);
        dum = ylim;
        dum1=[dum(1):0.01:dum(2)];
        for ii=1:(dum4-1);
            if (dd(ii+1,2) ~= dd(ii,2))
                dum3 = ones(size(dum1))*(ii+1);
                figure(fignum); plot(dum3,dum1,'k');
            end
        end
    case 1
        figure(fignum);
        plot(dd(1:dum4,7:last),'b');
        hold on;
        dum = ylim;
        dum1=[dum(1):0.01:dum(2)];
        for ii=1:(dum4-1);
           if (dd(ii+1,2) ~= dd(ii,2))
               dum3 = ones(size(dum1))*(ii+1);
               figure(fignum); plot(dum3,dum1,'k');
           end  
        end
        figure(fignum+1);
        plot(dd(dum4+1:last2,7:last),'r');
        hold on;
        for ii=dum4+1:(last2-1);
        if (dd(ii+1,2) ~= dd(ii,2))
               dum3 = ones(size(dum1))*(ii+1-dum4);
               figure(fignum+1); plot(dum3,dum1,'k');
           end  
        end
end






hold off;


