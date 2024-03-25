clc ;clear ;close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 代码来源：外文学位论文 Modelling and Coverage Improvement of DVB-T Networks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT PARAMETERS
% First, the required data for the program΄s operation are imported. These are the 
% following :
% a) Frequency f of the transmitted signal in MHz
% b) Transmitted Power (EIRP) in dBm.
% c) Transmitter Antenna height, Ht in meters
% d) Receiver antenna height, Hr in meters
% e) Transmitter Name
% f) Transmitter Longitude
% g) Transmitter Latitude
% h) Receiver name
% i) Transmitter Longitude
% j) Transmitter Latitude

% disp('INPUT PARAMETERS');
% disp('k-factor =1.33 (effect of the atmosphere) for all calculations.');
% % the above disp command works only for underlined blue text
% f=input('Operating Frequency f in MHz, f= ');
% Pt=input('EIRP in dBm, Pt= ');
% Ht=input('Antenna Transmitter height Ht in m, Ht= ');
% Hr=input('Antenna Receiver height Hr in m, Hr= ');
% Tname = input('Transmitter name? :','s');%'s' indicates that Tname is a string.
% % Transmitter coordinates
% x(1)=input('Transmitter LONG, LONG= ');
% y(1)=input('Transmitter LAT, LAT= ');
% Rname = input('Receiver name? :','s');%'s' indicates that Rname is a string.
% % Receiver coordinates
% x(2)=input('Receiver LONG, LONG= ');
% y(2)=input('Receiver LAT, LAT= ');

% 论文104页的
f = 217.25;
Pt = 79.18;
Ht = 45;
Hr = 2.5;
Tname = 'fashe';%'s' indicates that Tname is a string.
% Transmitter coordinates
x(1)=22.42255;
y(1)=41.80293;
Rname = 'jieshou';%'s' indicates that Rname is a string.
% Receiver coordinates
x(2)=22.488190;
y(2)=41.312160;



%% READ OF readhgt AND CREATING GEOPHYSICAL MAP AND ELEVATION PROFILE
% 阅读并创建地球物理图和高程剖面图
X=readhgt(39:41,21:24,'merge','crop','plot');
xi=linspace(x(1),x(2),4000);
yi=linspace(y(1),y(2),4000);
zi=interp2(X.lon,X.lat,double(X.z),xi,yi,'linear');
% Plot the red line in geophysical map
hold on,plot(xi,yi,'r'),hold off
% Display the names of transmitter and receiver on map in red line
txt1 = Tname;
text(x(1),y(1),txt1)
txt2 = Rname;
text(x(2),y(2),txt2)

%% 大圆的等矩形近似（Haversine公式）（RHUMB线或LOXODROME）
% Distance along this path, see equirectangular approximation.doc.
% I changed 6370 km erath radius in 6370000 in meters so the result is in
% meters,because MATLAB's default is in meters.
di=sqrt(((xi-x(1)).*cosd((y(1)+yi)/2)).^2+(yi-y(1)).^2)*6370000*pi/180;
% Beauducel gives the approximate formula below which gives similar results
% di=sqrt(((xi-x(1)).*cosd(yi)).^2+(yi-y(1)).^2)*6370000*pi/180;
% if you give here plot(di,zi)plots elevation profile without earth's bulge
dmin=di(1);
dmax=di(4000);
zmin=zi(1)+Ht;% Ht is trasmitter antenna height
Htransmitter=['Total Transmitter Antenna Height ',num2str(zmin),'m'];
disp(Htransmitter)
zmax=zi(4000)+Hr;% Hr is receiver antenna height.
Hreceiver=['Total Receiver Antenna Height ',num2str(zmax),'m'];
disp(Hreceiver)

%% 地球曲率（地球凸起）和大气衍射的校正
% In fact the effective curvature of the Earth can be calculated by adding
% a parabolic approximation of the effective Earth curvature to the path
% profile terrain data before finding horizons. One such parabolic function
% that can be used is r=dt*dr/2r=dt*(d-dt)/2r=dr*(d-dr)/2r. Where r=k*ro 
% and ro=earth radius=6370km=6370000m, and usually k=4/3=1.333.
% In MATLAB all parameters must be in meters(m),(default in MATLAB):
% so r=(di-dmin)*(dmax-di)/2*1.333*6.370.000=0.000000059*(di-dmin)*(dmax-di)
% See for further details, folder "Earth bulge" in PhD folder.
% Here : (di-dmin) is the distance between correction point and transmitter in meters
% and (dmax-di) is the distance between correction point and receiver in
% meters. The elevation changes if you change the constant k.
% So finally my correction formula is : 0.000000059*(di-dmin)*(dmax-di)
r=0.000000059*(di-dmin).*(dmax-di);
zk=zi+r; % adding Earth's bulge

%% 绘制立面纵断面
figure
hold on,plot(di,zk)
% Display the names of transmitter and receiver on elevation profile
% we multilpy Ht and Hr by 2, 2*Ht and 2*Hr to move text up, othewise it
% conflicts with the lines.
txt1 = Tname;
text(di(1),zi(1)+2*Ht,txt1)
txt2 = Rname;
text(di(4000),zi(4000)+2*Hr,txt2)
title(['Elevation Profile from ' sprintf('%s',char(Tname)) ' to ' sprintf('%s',char(Rname))]);
xlabel('distance in meters between stations');
ylabel('altitude in meters');
% plot the Earth curvature curve
hold on,plot(di,r,'g')
% plot the LOS red line in elevation profile
hold on,plot([dmin dmax],[zmin zmax],'r')
% plot the antennas in figure
% Transmitter Antenna Height, the formula is: line([X X], [LowY HighY])
% We plot again Ht to show Transmitter Antenna, although it was plotted in zmin=zi(1)+Ht
line([di(1) di(1)],[zi(1) zi(1)+Ht],'color','k','LineWidth',2)
% Receiver Antenna Height, the formula is: line([X X], [LowY HighY])
% We plot again Hr to show Receiverer Antenna, thow it was plotted in
% zmax=zi(4000)+Hr
line([di(4000) di(4000)],[zi(4000) zi(4000)+Hr],'color','k','LineWidth',2)

%% 程序绘制了第一菲涅耳区和0.6倍第一菲涅耳区，并与0.6倍第一菲涅尔区的高程剖面相交。
x1=dmin;
y1=zmin;
x2=dmax;
y2=zmax;
% Plot of the 1st Fresnel zone
fr=f*1000000;% f in Hz218
c=300000000;% c light speed 300000000m/s
l=c/fr; % l=wave length in meters
a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
r = sqrt(l*a/2);% b=r
t = linspace(0,2*pi);
X = a*cos(t);
Y = r*sin(t); 
w = atan2(y2-y1,x2-x1);
x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
hold on, plot(x,y,'y')
grid on
% Plot the 0.6 times Fresnel zone, the upper half which does not intersects
% with the obstacle
r = 0.6*sqrt(l*a/2);% b=r
t = linspace(0,2*(pi/2));
X = a*cos(t);
Y = r*sin(t); 
w = atan2(y2-y1,x2-x1);
x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
hold on, plot(x,y,'m')
% Plot the 0.6 times Fresnel zone, the bottom half which intersects
% with the obstacles
r = 0.6*sqrt(l*a/2);% b=r
t = linspace(0,2*(-pi/2));
X = a*cos(t);
Y = r*sin(t); 
w = atan2(y2-y1,x2-x1);
x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
hold on, plot(x,y,'m')
% INTERSECT THE 0.6F AND THE ELEVATION PROFILE.
[ki,li]=polyxpoly(x, y, di, zk, 'unique');% polyxpoly from mapping Toolbox.
% 'unique' filters out duplicate intersections, which may result 
% if the input polylines are self-intersecting.
mapshow(ki, li, 'DisplayType','point','Marker','o');

%% 程序检查是否存在交叉点，如果没有交叉点，则计算自由空间路径损耗和场强，
%  如果存在两个交叉点（意味着一个障碍物），则程序计算位于该点的接收器处的单刀
%  边缘路径损耗和光强。
%   CHECK IF MATRIX [ki,li] IS EMPTY (FREE SPACE PATH LOSS) OR NOT.
%====================================================================
if isempty([ki,li]) % FREE SPACE PATH LOSS CALCULATION
disp('No obstructions detected')
% COMPUTE FREE SPACE PATH LOSS IN DBI WHERE f in MHZ and d in km
FSPL=32.45+20*log10(f)+20*log10(di(4000)/1000);
V=['Free Space Path Loss=',num2str(FSPL),'dB'];
disp(V)
% CALCULATION: E(dBμV/m)=P(dBm)-PL+20logf(MHz)-Gr(dBi)+77.2
% IN OUR COMPUTE WE CONSIDER Gr(dBi)=0
E=Pt-FSPL+20*log10(f)+77.2;
U=['Field Strength E =',num2str(E),' dBuV/m'];
disp(U)
%====================================================================
else
% Calculation of rows and columns of matrix [ki,li]
% sizeOfMatrix = size([ki,li])
% Calculation of rows of matrix [ki,li], ki-->Distance, li-->Height
nRows = size([ki,li], 1); % 1 stands for the first dimension
Row=nRows/2;% because each obstructions is counted twice, we devide by 2
disp(' Distance(ki) Height(li) (in meters) ')
disp([ki,li])
M=['There are ',num2str(Row),' obstructions at ki and li positions.'];
disp(M)
DASH='===============================================================';
disp(DASH)
DASH='===============================================================';
disp(DASH)
%====================================================================
if Row==1
%====================================================================
% If Row=1 MEANS THAT WE HAVE 1 OBSTACLE SO SINGLE KNIFE-EDGE MODEL
%====================================================================
disp('Single Knife-Edge calculation begins');
% I consider that I have only 2 intersection points li(1)and li(2)
% because in Knife-Edge model we take into account the first obstacle
% otherwise the model becomes too complicated.li-->Heihgt, ki-->distance
% Find the indices of c1 and e1 in matrix di.
index=find(di<ki(2) & di>ki(1));
% Find the first and last element of indices matrix index, 
% index(1) and index(end).
f1=index(1);
f2=index(end);
% Using indices from above find the elements in zk creating a submatrix mi
mi=zk(f1:f2);
% find the max of matrix mi, maxheight.
maxheight=max(mi);
Mh=['Height of Single Knife-Edge=',num2str(maxheight),' m'];
disp(Mh)
% find the index of maxheight in matrix mi which will 
% be the same in matrix d1,sometimes we have two (2) same index 
% so we take the first one, index(1).
% In general, when we are dealing with floats, 
% always try to avoid exact comparisons and look for proximity,
% that's why we use the formula below instead of index=find(mi==maxheight).
index = find( abs(mi-maxheight) < 0.0001 );
f3=index(1);
% Create a submatrix d1 of di with indices f1,f2.
d1=di(f1:f2);
D1=d1(f3);
Md=['Horizontal distance between Transmitter and Single Knife-Edge =',num2str(D1),' m'];
disp(Md)
% Now I find Earth radius r1 in this point because I need the coordinates of
% the line "maxheight" which are : start point(D1,r1)and
% end point(D1,maxheight) to plot Knife-Edge1 with the use of
% the command "line"
r1=0.000000059*(D1-dmin).*(dmax-D1);
% Plot the vertical line of Single Knife-Edge
line([D1 D1],[r1 maxheight],'color','k','LineWidth',1)
% Intersect the 1st Knife-Edge line and the LOS Line
% The command "polyxpoly" does not work here because these two lines
% don't intersect but the their extensions intersect, so if I use
% the command "polyxpoly" i get the message 'empty matrix',that's why we 
% use algebra in finding intersection points between two lines,221
% see explanation doc.
% In MATLAB backslash operator (\)is used to solve systems of linear
% equations Ax=B ==>x=A\B.
A=[(zmax-zmin),-(dmax-dmin);(maxheight-r1),-(D1-D1)];
b=[(zmax-zmin)*dmin-(dmax-dmin)*zmin;(maxheight-r1)*D1-(D1-D1)*r1];
Intersectionpoint=A\b;
% Parameter h in formula v, i.e. height between Knife-Edge and LOS line
h=maxheight-Intersectionpoint(2);
Ph=['Height between Single Knife-Edge line and LOS, h =',num2str(h),' m'];
disp(Ph)
%---------------------------------------------------------------------------------------------
% Find the length of the LOS red line.
length = sqrt((dmax-dmin).^2+(zmin-zmax).^2);
LOS=['LOS length =',num2str(length),' m'];
disp(LOS)
%---------------------------------------------------------------------------------------------------------
% Find the distance between Knife-Edge and Receiver
% which is D2=dmax-D1.
D2=dmax-D1;
DHrKE=['Horizontal distance between Single Knife-Edge and Receiver =',num2str(D2),' m'];
disp(DHrKE)
% Calculating v parameter of fresnel Integral
v=h*sqrt((2*dmax)/(l*D1*D2));
Fv=['Fresnel Parameter v=',num2str(v)];
disp(Fv)
% Approximate solution for equation Gd(dB)=20log[F(v)]provided by Lee.
Gd=KnifelossLee(v);
KED1=['Single Knife-Edge Path Loss Gd =',num2str(Gd),' dB'];
disp(KED1)
%---------------------------------------------------------------------------------------------------------
% Computing the Total Path Loss after Single Knife-Edge Model
FSPL=32.45+20*log10(f)+20*log10(di(4000)/1000);
V=['Free Space Path Loss =',num2str(FSPL),' dB'];
disp(V)
% The Total Path Loss in dB is the sum of Free Space Path Loss 
% and Knife-Edge Path Loss.
Gtotal=FSPL+abs(Gd);
Gt=['Single Knife-Edge Path Loss =',num2str(Gtotal),' dB'];
disp(Gt)
% And in dBuV/m is :
Etotal=Pt-Gtotal+20*log10(f)+77.2;
Et=['Field Strength =',num2str(Etotal),' dBuV/m'];
disp(Et)
% END OF SINGLE KNIFE EDGE CALCULATION
 else
% AUTOMATION PROGRAM FINDS ALL THE HEIGHTS AND MAX1 AND MAX2 FOR 
% FRESNEL PARAMETERS V1 & V2
% Preallocation of matrix V
V=zeros(1,nRows);
%-----------------------------------------------------------------------------------------------------------------
for a=1:2:nRows
indexA=find(di<ki(a+1) & di>ki(a));
% Find the first and last element of indices matrix index, 
% index(1) and index(end).
f1=indexA(1);
f2=indexA(end);
% Using indices from above find the elements in zk creating a submatrix mi
mi=zk(f1:f2);
% find the max of matrix mi, maxheight.
maxheight=max(mi);
Mh=['Height of Knife-Edge=',num2str(maxheight),' m'];
disp(Mh)
% find the index of maxheight in matrix mi which will 
% be the same in matrix d1,sometimes we have two (2) same index 
% so we take the first one,index(1).
% In general, when we are dealing with floats, 
% always try to avoid exact comparisons and look for proximity,
% that's why we use the formula below instead of index=find(mi==maxheight).
indexB = find( abs(mi-maxheight) < 0.0001 );
f3=indexB(1);
% Create a submatrix d1 of di with indices f1,f2.223
d1=di(f1:f2);
distance=d1(f3);% distance is D1
Md=['Horizontal distance of Knife-Edge from transmitter =',num2str(distance),' m'];
disp(Md)
% Earth radius r1 in this point.
r1=0.000000059*(distance-dmin).*(dmax-distance);
line([distance distance],[r1 maxheight],'color','k','LineWidth',1)
% Intersect the Knife-Edge line and the LOS Line
A=[(zmax-zmin),-(dmax-dmin);(maxheight-r1),-(distance-distance)];
b=[(zmax-zmin)*dmin-(dmax-dmin)*zmin;(maxheight-r1)*distance-(distance-distance)*r1];
Intersectionpoint=A\b;
% Parameter h in formula v, i.e. height between Knife-Edge line and LOS
h=maxheight-Intersectionpoint(2);
Ph=['Height between Knife-Edge line and LOS, h =',num2str(h),' m'];
disp(Ph)
% Find the length of the LOS red line.
length = sqrt((dmax-dmin).^2+(zmin-zmax).^2);
LOS=['LOS length =',num2str(length),' m'];
disp(LOS)
% Find distance D1 between Transmitter and Knife-Edge.
D1=distance;
DHtKE=['Horizontal distance between Transmitter and Knife-Edge =',num2str(D1),' m'];
disp(DHtKE)
% Find distance D2 between Knife-Edge and Receiver.
D2=dmax-D1;
DHrKE=['Horizontal Distance between Knife-Edge and Receiver =',num2str(D2),' m'];
disp(DHrKE)
% Calculating v parameter of fresnel Integral
v=h*sqrt((2*(D1+D2))/(l*D1*D2));
Fv=['Fresnel Parameter v=',num2str(v)];
disp(Fv)
% Here we create a matrix V(a), that contains all the parameters v from the
% for...next loop
V(a)=h*sqrt((2*(D1+D2))/(l*D1*D2));
DASH='==============================================================';
disp(DASH)
% Because the matrix V(a) contains zeros and so it gives wrong indexes of
% the maximum value 1 (max1) and the maximum value 2 (max2), we delete all
% zeros from matrix V(a) creating a matrix L that has all the data of
% matrix V(a) but without zeros.224
L = V(V~=0);
[max1,idx1]= max(L);% we find the max value.
L(idx1)=NaN;
[max2,idx2]= max(L);% we find the second max value.
L(idx1)=max1;
DASH='==============================================================';
disp(DASH)
end
DASH='=============================================================';
disp(DASH)
MAXV=['The maximum V1 =',num2str(max1),' and it is caused by obstacle ' ,num2str(idx1)];
disp(MAXV)
MAXV=['The maximum V2 =',num2str(max2),' and it is caused by obstacle ' ,num2str(idx2)];
disp(MAXV)
DASH='==============================================================';
disp(DASH)

%% 
%====================================================================
% DEYGOUT MODEL
%====================================================================
% 1st KNIFE-EDGE
%====================================================================
% li-->Heihgt, ki-->distance
% Find the indices of c1 and e1 in matrix di.
% The program gives the number of the obstacle, i.e 9 obstacle but we want
% the intersection points. For the 9th(idx1=9) obstacle the intersection 
% points are 17 and 18.How do we find these? We say a2=2*idx1=2*9=18 and
% a1=a2-1=18-1=17. That is the logic for lines 345 ,346 (below),422 and 423.
a2=2*idx1;
a1=a2-1;
index=find(di<ki(a2) & di>ki(a1));
% Find the first and last element of indices matrix index, 
% index(1) and index(end).
f1=index(1); 
f2=index(end);
% Using indices from above find the elements in zk creating a submatrix mi
mi=zk(f1:f2);
% Find the max of matrix mi.
maxheight1=max(mi);
Mh=['Height of 1st Knife-Edge=',num2str(maxheight1),' m'];
disp(Mh)
% find the index of maxheight1 in matrix mi which will 
% be the same in matrix d1,sometimes we have two (2) same index 
% so we take the first one,index(1).
% In general, when we are dealing with floats, 
% always try to avoid exact comparisons and look for proximity,
% that's why we use the formula below instead of index=find(mi==maxheight1).
index = find( abs(mi-maxheight1) < 0.0001 );
f3=index(1);
% Create a submatrix d1 of di with indices f1,f2.
d1=di(f1:f2);
D1=d1(f3);
%--------------------------------------------------------------------------
% Earth radius r1 in this point.
r1=0.000000059*(D1-dmin).*(dmax-D1);
%--------------------------------------------------------------------------
% Plot the vertical line of the 1st Knife-Edge.
line([D1 D1],[r1 maxheight1],'color','k','LineWidth',1)
% Intersect the Knife-Edge line and the LOS Line
A=[(zmax-zmin),-(dmax-dmin);(maxheight1-r1),-(D1-D1)];
b=[(zmax-zmin)*dmin-(dmax-dmin)*zmin;(maxheight1-r1)*D1-(D1-D1)*r1];
Intersectionpoint=A\b;
% Parameter h1 in formula v, i.e. height between 1st Knife-Edge and LOS
h1eff=maxheight1-Intersectionpoint(2);
Ph=['Height between 1st Knife-Edge and LOS, h1eff =',num2str(h1eff),' m'];
disp(Ph)
% Find the length of the LOS red line.
length = sqrt((dmax-dmin).^2+(zmin-zmax).^2);
LOS=['LOS length =',num2str(length),' m'];
disp(LOS)
% Distance between Transmitter and 1st Knife-Edge,
DHtKE=['Horizontal distance between Transmitter and 1st Knife-Edge D1 =',num2str(D1),' m'];
disp(DHtKE)
% Distance between 1st Knife-Edge and Receiver.
D2=dmax-D1;
DHrKE=['Horizontal distance between 1st Knife-Edge and Receiver D2=',num2str(D2),' m'];
disp(DHrKE)
% Calculating v1 parameter of fresnel Integral
v=h1eff*sqrt((2*dmax)/(l*D1*D2));
Fv=['Fresnel Parameter v1=',num2str(v)];
disp(Fv)
% Approximate solution for equation Gd(dB)=20log[F(v)]provided by Lee.
Gd1=KnifelossLee(v);
SKE=['1st Knife-Edge Path Loss Gd1 =',num2str(Gd1),' dB'];
disp(SKE)
DASH='---------------------------------------------------------------------';
disp(DASH)
%====================================================================
% 2nd KNIFE-EDGE
%====================================================================
% 2nd Knife Edge calculation begins
% li-->Heihgt, ki-->distance
% Find the indices of c1 and e1 in matrix di.
a4=2*idx2;
a3=a4-1;
index=find(di<ki(a4) & di>ki(a3));
% Find the first and last element of indices matrix index, 
% index(1) and index(end).
f1=index(1);
f2=index(end);
% Using indices from above find the elements in zk creating a submatrix mi
mi=zk(f1:f2);
% find the max of matrix mi, maxheight2.
maxheight2=max(mi);
Mh=['Height of 2nd Knife-Edge=',num2str(maxheight2),' m'];
disp(Mh)
% find the index of maxheight2 in matrix mi which will 
% be the same in matrix d1,sometimes we have two (2) same index 
% so we take the first one,index(1).
% In general, when we are dealing with floats, 
% always try to avoid exact comparisons and look for proximity,
% that's why we use the formula below instead of index=find(mi==maxheight2).
index = find( abs(mi-maxheight2) < 0.0001 );
f3=index(1);
% Create a submatrix d1 of di with indices f1,f2.
d1=di(f1:f2);
D3=d1(f3);
Md=['Distance between Transmitter and 2nd Knife-Edge D3=',num2str(D3),' m'];
disp(Md)
% Find distance between 2nd Knife-Edge and Receiver.
D4=abs(dmax-D3);
DHrKE2=['Distance between Receiver and 2nd Knife-Edge D4=',num2str(D4),' m'];
disp(DHrKE2)
% Distance between the two Knife-Edge is D5.
D5=abs(D3-D1);
DM1M2=['Distance between the Two Knife-Edges D5=',num2str(D5),' m'];
disp(DM1M2)
%--------------------------------------------------------------------------
% Earth radius r2 in this point.
r2=0.000000059*(D3-dmin).*(dmax-D3);
%--------------------------------------------------------------------------
% Plot the vertical line of the second Knife-Edge.
line([D3 D3],[r2 maxheight2],'color','k','LineWidth',1)
% Intersect the Knife-Edge 2 line and the LOS Line.
A=[(zmax-zmin),-(dmax-dmin);(maxheight2-r2),-(D3-D3)];
b=[(zmax-zmin)*dmin-(dmax-dmin)*zmin;(maxheight2-r2)*D3-(D3-D3)*r2];
Intersectionpoint=A\b;
% Parameter h in formula v, i.e. height between Knife-Edge line and LOS
h2=maxheight2-Intersectionpoint(2);
Ph2=['Height between 2nd Knife-Edge and LOS, h2 =',num2str(h2),' m'];
disp(Ph2)
%--------------------------------------------------------------------------
% Plot the line from maxheight1 to Transmitter, i.e. line TM1
% Coordinates are (D1,maxheight1) and (dmin,zmin).
line([dmin D1],[zmin maxheight1],'color','k')
%--------------------------------------------------------------------------
% Plot the line from maxheight1 to Receiver, i.e. line M1R.
% Coordinates are (D1,maxheight1) and (dmax,zmax).
line([dmax D1],[zmax maxheight1],'color','k')
%--------------------------------------------------------------------------
% HERE WE CALCULATE INTERSECTION POINT BETWEEN MAXHEIGHT2 AND M1R 
% IF SECONDARY OBSTACLE LIES ON THE RIGHT OF PRIMARY OR M1T IF 
% SECONDARY OBSTACLE LIES ON THE LEFT OF PRIMARY
% Observe that h2eff is the same with Epstein-Peterson h2 and Giovaneli h2.
%=====================================================================
% IF D1<D3 MEANS THAT THE SECONDARY OBSTACLE LIES ON THE RIGHT OF THE MAIN
%====================================================================
if D1 < D3
% D1<D3 means that secondary obstacle M2 lies on the right of the main
% obstacle M1, so its maxheight2 intersects with line M1R.
% So we intersect the line maxheight1-R(M1R) and maxheight2
% Coordinates are (distance1,maxheight1) and (distance2,maxheight2).
% we use distances D5,D4 for the secondary obstacle.
A=[(maxheight1-zmax),-(D1-dmax);(maxheight2-r2),-(D3-D3)];
b=[(maxheight1-zmax)*dmax-(D1-dmax)*zmax;(maxheight2-r2)*D3-(D3-D3)*r2];
Intersectionpoint=A\b;
h2eff=maxheight2-Intersectionpoint(2);
DH2=['Height between maxheight2 and line M1R h2eff=',num2str(h2eff),' m'];
disp(DH2)
D5=abs(D3-D1);
D4=abs(dmax-D3);
% calculating v2
v=h2eff*sqrt((2*(D4+D5))/(l*D4*D5));
Fv=['Fresnel Parameter v2=',num2str(v)];
disp(Fv)
% Approximate solution for equation Gd(dB)=20log[F(v)]provided by Lee.
Gd2=KnifelossLee(v);
DKE=['2nd Knife-Edge Path Loss Gd2 =',num2str(Gd2),' dB'];
disp(DKE)
DASH='---------------------------------------------------------------------';
disp(DASH)
%====================================================================
% ELSE D1>D3 MEANS THAT THE SECONDARY OBSTACLE LIES ON THE LEFT OF THE MAIN
%====================================================================
else
% else means that D1>D3 and the secondary obstacle M2 lies on the left of the main
% obstacle, so its maxheight2 intersects with line M1T.
% So we intersect the line maxheight1-T(M1T) and maxheight2.
% Coordinates are (distance1,maxheight1) and (distance2,maxheight2).
% Now we use distances D3,D5.
A=[(maxheight1-zmin),-(D1-dmin);(maxheight2-r2),-(D3-D3)];
b=[(maxheight1-zmin)*dmin-(D1-dmin)*zmin;(maxheight2-r2)*D3-(D3-D3)*r2];
Intersectionpoint=A\b;
h2eff=maxheight2-Intersectionpoint(2);
DH2=['Height between maxheight2 and line M1T h2eff=',num2str(h2eff),' m'];
disp(DH2)
% calculating v2229
v=h2eff*sqrt((2*(D3+D5))/(l*D3*D5));
Fv=['Fresnel Parameter v2=',num2str(v)];
disp(Fv)
% Approximate solution for equation Gd(dB)=20log[F(v)]provided by Lee.
Gd2=KnifelossLee(v);
DKE=['2nd Knife-Edge Path Loss Gd2 =',num2str(Gd2),' dB'];
disp(DKE)
DASH='---------------------------------------------------------------------';
disp(DASH)
end
%====================================================================
% TOTAL PATH LOSS AFTER TWO KINFE-EDGES
%====================================================================
FSPL=32.45+20*log10(f)+20*log10(di(4000)/1000);
V=['Free Space Path Loss =',num2str(FSPL),' dB'];
disp(V)
% The Total Path Loss in dB is the sum of Free Space Path Loss 
% and the 2 Knife-Edge Path Loss, Gd1 and Gd2.
Gtotal=FSPL+abs(Gd1)+abs(Gd2);
Gt=['Double Knife-Edge Path Loss (DEYGOUT) =',num2str(Gtotal),' dB'];
disp(Gt)
% And in dBuV/m is :
Etotal=Pt-Gtotal+20*log10(f)+77.2;
Et=['Field Strength =',num2str(Etotal),' dBuV/m'];
disp(Et)
DASH='*********************************************************************';
disp(DASH)
MINFSA=[' The minimum Field Strength is ',num2str(round(Etotal,1)),' dBμV/m .'];
disp(MINFSA)
MINFSB=['and is caused between the Main Obstacle ',num2str(idx1),' and the Secondary obstacle ',num2str(idx2), '.',];
disp(MINFSB)
DASH4='*********************************************************************';
disp(DASH4)
end
end 












