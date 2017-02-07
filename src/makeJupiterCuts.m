function [data] = makeJupiterCuts(filename,degreesAwayFromCenter)
tic
energyName = regexp(filename,'_','split');
energyName=energyName{3};

fid = fopen(filename);
mydata = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %*[^\n]', 'delimiter', ',','CollectOutput',1);
data=mydata{1};
data=[data(:,1:5) data(:,10)];
fclose(fid);

RA=data(:,2);
DEC=data(:,3);
photonTime=data(:,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Get RA/DEC functions %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('ephem_jup_2008-2014.txt');

mydata = textscan(fid, '%f %f %f %[^\n}]', 'delimiter', ',','CollectOutput',1);
jupiterData =mydata{1};
jupiterData(1,:)=[];           % Remove first line
fclose(fid);
jupiterData(:,2)=(360/24).*jupiterData(:,2);
jupiterTime=jupiterData(:,1);
jupiterTime = (jupiterTime - 51910).*86400;
position=jupiterData(:,2:3);

RAfunction = fit(jupiterTime,position(:,1),'fourier8');
DECfunction = fit(jupiterTime,position(:,2),'fourier8');
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now create the cuts!
% For each source, go through all the photons

[rows,cols]=size(data);

jupiterRow=1;
nonJupiterRows=1;
jupiterPhotons=zeros(1,6);
jupiterPosition=zeros(rows,2);
for i=1:rows
    % For each photon, plug time into functions and see if < _ degrees away
    current_time=photonTime(i);
    jupiterPosition(i,:)=[RAfunction(current_time) DECfunction(current_time)];
    
    if angDist([RA(i) DEC(i)],jupiterPosition(i,:)) < degreesAwayFromCenter
        jupiterPhotons(jupiterRow,:)=data(i,:);
        jupiterRow=jupiterRow+1;
    else
        nonJupiterPhotonRows(nonJupiterRows,1)=i;
        nonJupiterRows=nonJupiterRows+1;
    end
        
end

loc = ['jupiterData_' energyName '_' num2str(degreesAwayFromCenter) '_deg.mat'];
save(loc,'jupiterPhotons','nonJupiterPhotonRows','jupiterPosition');

fid = fopen(['Non-jupiter_rows_' energyName '_' num2str(degreesAwayFromCenter) '_deg.txt'], 'wt');
fprintf(fid,'%i\n',nonJupiterPhotonRows);
fclose(fid);
toc
end

