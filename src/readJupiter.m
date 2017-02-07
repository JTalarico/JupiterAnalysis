clear all

fid = fopen('ephem_jup_2008-2014.txt');

mydata = textscan(fid, '%f %f %f %[^\n}]', 'delimiter', ',','CollectOutput',1);
data =mydata{1};
data(1,:)=[];           % Remove first line
fclose(fid);
data(:,2)=(360/24).*data(:,2);
time=data(:,1);

time = (time - 51910).*86400;

position=data(:,2:3);

[rows,~]=size(position);
foo=zeros(rows,1);
for i=1:rows
    foo(i,1) = time(i)-time(1);
    foo(i,2) = angDist(position(1,:),position(i,:));
end

RAfunction = fit(time,position(:,1),'fourier8');
DECfunction = fit(time,position(:,2),'fourier8');
           
           
           
test(:,1)= RAfunction(time)-position(:,1);
test(:,2)=DECfunction(time)-position(:,2);

%DECfunction = @(Xtime) ().*Xtime.^8 + ().*Xtime.^7 + ().*Xtime.^6 + ().*Xtime.^5 + ().*Xtime.^4 + ().*Xtime.^3 +().*Xtime.^2 +()*Xtime + ();