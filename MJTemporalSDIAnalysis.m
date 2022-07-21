%%%%%%%%%%%%%%%%%%% SDI Count Work%%%%%%%%%%%%%%%%%%%%%%%%

%% This work reshuffles the success tracked particles
%Get the number of iterations to run from the app.track

MJnumcolumn_m2 = size(app.init_track);


%Result usually comes through as a double, we just want the end value (max)
MJnumcolumn_m = max(MJnumcolumn_m2);

%Create a blank table ready to store results to
MJtable = cell(1);


%% Run through each iteration and combine the start x,y coordinate points
%with the end x,y points, as well as its success (0,1, or 2)
%Results are stored as a cell named (MJcell)
for i = 1 : MJnumcolumn_m
    
    MJa = (app.init_track{i});
    MJb = (app.fin_track{i});
    MJc = (app.success_track{i});
    
    
    MJd = [MJa, MJb, MJc];   
    MJcell{i} = {MJd};
    
end

%% When plotting as a scatter, it helps to have the background points
%included to help with visualisation. This just takes the last run of the
%iterations and creates its own variable ready to pull the background
%points from.
background = MJd;
background(any((background(:,5)==1),2), : ) = [];
background(any((background==0),2), : ) = [];

%% Turn the cell MJcell into a table ready for reading. Run this as a command
%to reset the MJtable
MJtablesuccess = cell2table(MJcell);
MJtablefailed = cell2table(MJcell);

%% Go through the table and remove the values associated with the "==" value
%seen below (0,1, or 2), depending on what values you want to leave. This
%is for successfull tracers
for i = 1 : MJnumcolumn_m
    MJtablesuccess{1,i}{1,1}(any((MJtablesuccess{1,i}{1,1}==0),2), : ) = [];
    MJtablesuccess{1,i}{1,1}(any((MJtablesuccess{1,i}{1,1}==2),2), : ) = [];
end

%% Go through the table and remove the values associated with the "==" value
%seen below (0,1, or 2), depending on what values you want to leave. This
%is for unsuccessfull tracers
for i = 1 : MJnumcolumn_m
    MJtablefailed{1,i}{1,1}(any((MJtablefailed{1,i}{1,1}==1),2), : ) = [];
    MJtablefailed{1,i}{1,1}(any((MJtablefailed{1,i}{1,1}==2),2), : ) = [];
end

%% NOTE : Outputs

%The outputs are MJtablesuccess that has all of the initial positions of
%tracers that have been successfully tracked.


%%
%%%%%%%%%%%%%%%%%%% SDI Region Work%%%%%%%%%%%%%%%%%%%%%%%%

% Define a uniform ROI within KLT's ROI

xy = app.boundaryLimitsM(:,1:2);
x = xy(:,1);
y = xy(:,2);

%% Temporal Analysis
outside = 1;  
grid_val = app.ResolutionmpxEditField.Value*app.BlocksizepxEditField.Value  ;

% Mesh grid the axis into the grid value defined and draw on the ROI points
% given in KLT
[X,Y] = meshgrid(min(x)-outside:grid_val:max(x)+outside, min(y)-outside:grid_val:max(y)+outside);    

%% Deifne the ROI, Has to start with the top left for this to work 
% and go clockwise
%y1 Bottom Right to bottom left
m1=(y(4)-y(1))/(x(4)-x(1));
%y2 top right to bottom right
m2=(y(1)-y(2))/(x(1)-x(2));
%y3 Top left to top right
m3=(y(2)-y(3))/(x(2)-x(3));
%y4 bottom left to top left
m4=(y(3)-y(4))/(x(3)-x(4));

%Bottom
Ytest1=m1*(X(1,:)-x(1))+y(1);

%RHS
Ytest2=m2*(X(1,:)-x(2))+y(2);

%Top
Ytest3=m3*(X(1,:)-x(3))+y(3);

%LHS
Ytest4=m4*(X(1,:)-x(4))+y(4);

%Comment these out one by one if you cant get it to work
YY=Y;
%Bottom
YY(YY<Ytest1)=NaN;
%RHS
YY(YY<Ytest2)=NaN;
%Top
YY(YY>Ytest3)=NaN;
%LHS
YY(YY>Ytest4)=NaN;

YY(~isnan(YY))=1;
mask=YY;

XXX=X.*mask;
YYY=Y.*mask;


%% Analyse the data to remove incomplete grid points and make them the 
%centre points of the analysis 

zcx = zeros(size(YY,1),1); %Create zero column
zcy = zeros(1,size(YY,2)); %Create zero row

%Shift to right by 1 grid
newmatrix = [zcx, YY]; %add to beginning
newmatrix(:,size(newmatrix,2)) = []; %delete last column 
%Delete outside points
true1 = newmatrix==YY;
true1 = +true1;
true1(true1==0) = NaN;

%Shift to left
newmatrix = [true1,zcx]; %add to end
newmatrix(:,1) = []; %delete first column 
%Delete outside points
true2 = newmatrix==true1;
true2 = +true2;
true2(true2==0) = NaN;

%Shift up
newmatrix = [zcy; true2]; %add to top
newmatrix(size(newmatrix,1),:) = []; %delete last row 
%Delete outside points
true3 = newmatrix==true2;
true3 = +true3;
true3(true3==0) = NaN;

%Shift down
newmatrix = [true3; zcy]; %add to bottom
newmatrix(1,:) = []; %delete first row 
%Delete outside points
true4 = newmatrix==true3;
true4 = +true4;
true4(true4==0) = NaN;

newXXX = X.*true4;
newYYY = Y.*true4;

regionofinterest = figure('Name', 'Region of Interest','Color','w')
[X,Y] = meshgrid(min(x)-outside:grid_val:max(x)+outside, min(y)-outside:grid_val:max(y)+outside);   
hold on 
hm = mesh(X,Y,X*0);  
hp = plot([x],[y],'r-') 
set(hm,'EdgeColor','k')  
set(hp,'LineWidth',2)  
set(gca,'Visible','off')
hold off
%     
figure(regionofinterest)
hold on
plot(XXX,YYY,'or','MarkerFaceColor','r','MarkerSize',5)
figure(regionofinterest)
hold on
plot(newXXX,newYYY,'or','MarkerFaceColor','g','MarkerSize',5)
hold off
% 
%% Need to start segregating the data into its grid
MeanGridDensity = nan(1,MJnumcolumn_m2(2));
TracePointDensity = MeanGridDensity;
CV_Density = MeanGridDensity;
Dispersion = MeanGridDensity;
NumParticles = []

for cont = MJnumcolumn_m2(1):MJnumcolumn_m2(2)
    W1 = MJtablesuccess{1,cont}{1}(:,1);
    W2 = MJtablesuccess{1,cont}{1}(:,2);
    
    
    
    V1 = rmmissing(newXXX(:));
    V2 = rmmissing(newYYY(:));
   
    histW = [W1,W2];
    N = hist3(histW,'Ctrs',{(min(V1):grid_val:max(V1)) min(V2):grid_val:max(V2)});
    
    true5 = true4 ;
    true5(isnan(true5)) = 0;
    true5( ~any(true5,2), : ) = [];  %rows
    true5( :, ~any(true5,1) ) = [];  %columns
    true5 = transpose(true5);
    N = true5.*N ;
    
    truex = transpose(newXXX);
    truey = transpose(newYYY);
    
    truex( ~any(truex,2), : ) = [];  %rows
    truex( :, ~any(truex,1) ) = [];  %columns
    truey( ~any(truey,2), : ) = [];  %rows
    truey( :, ~any(truey,1) ) = [];  %columns
    
    
    Xgrid = truex(:);
    Ygrid = truey(:);
%%
   %%%%%%%%%%%%%%%%%%%% Calculate Values Needed for SDI
    NumParticles = N(:);
    MeanGridDensity(cont) = (mean(N(:)))
    MinNumberInGrid = (nanmin(N(:)))
   
    DefinedResolution  = app.ResolutionmpxEditField.Value; %defined by user  
    
    AreaOfGrids = grid_val^2;
    NumberOfGrids = size(N(:));
    TotalArea = AreaOfGrids*(NumberOfGrids(1));
    PeriodicTracers = size(W1);
    ROIPixels = TotalArea/(DefinedResolution^2);  
    
    %IMPORTANT VALUE: DENSITY
    TracePointDensity(cont) = PeriodicTracers(1)/ROIPixels;
    %
    
    VarGridDensity = var(N(:))
    %COULD BE IMPORTANT?: Coeffecient of Variation
    CV_Density(cont) = sqrt(VarGridDensity)./MeanGridDensity(cont);
    %

    %IMPORTANT VALUE: DISPERSION
    Dispersion(cont)= VarGridDensity/MeanGridDensity(cont); % Computation of the Aggregation/Dispersion parameter D*
    %
    
end

MeanDensity = nanmean(TracePointDensity);
MeanCV_Area = nanmean(CV_Density);
MeanNu = nanmean(Dispersion);

% NEW SDI has been made, but only relates to individual videos to show the
% period of time that is ideal for processing. 
IdealDispersion = 1
MaxTracePointDensity = max(TracePointDensity)
IdealSDI = IdealDispersion./MaxTracePointDensity
MJSDI = ((Dispersion)./(TracePointDensity))./ IdealSDI;

%ALONSO'S SDI WORK:
SDI = ((Dispersion.^0.1)./(TracePointDensity./(1.02E-03)));

tensec = 10/ app.ExtractionratesEditField.Value(1);
MJSDImean = movmean(MJSDI,tensec ,'Endpoints','discard');
SDImean = movmean(SDI,tensec , 'Endpoints','discard');

[valMJ indMJ] = min(MJSDImean(:));
[valSDI indSDI] = min(SDImean(:));
[valMJbad indMJbad] = max(MJSDImean(:));
[valSDIbad indSDIbad] = max(SDImean(:));


%% Plotting
figure('Name', 'SDI Results')
subplot(5,1,1)
plot(MJnumcolumn_m2(1):(MJnumcolumn_m2(2)),TracePointDensity)
xlabel('Interval (Extraction Rate * s)')
ylabel('\rho (ppp)')

subplot(5,1,2)
plot(MJnumcolumn_m2(1):(MJnumcolumn_m2(2)),Dispersion)
xlabel('Interval (Extraction Rate * s)')
ylabel('D^*')

subplot(5,1,3)
plot(MJnumcolumn_m2(1):(MJnumcolumn_m2(2)),CV_Density)
xlabel('Interval (Extraction Rate * s)')
ylabel('CV Grid Density')

subplot(5,1,4)
plot(MJnumcolumn_m2(1):(MJnumcolumn_m2(2)),SDI)
hold on
scatter(indSDI,valSDI)
scatter(indSDIbad,valSDIbad)
hold off
xlabel('Interval (Extraction Rate * s)')
ylabel('SDI')

subplot(5,1,5)
plot(MJnumcolumn_m2(1):(MJnumcolumn_m2(2)),MJSDI)
hold on
scatter(indMJ,valMJ)
scatter(indMJbad,valMJbad)
hold off
xlabel('Interval (Extraction Rate * s)')
ylabel('MJSDI')


%%
%%%END OF ANALYSIS

% 
% %% Collates all of the x and y points of the successfull tracers into two
% %seperate variables for a full analysis on the whole video

WW1=[]
WW2=[]
if W1 == 0 
    for i = 1:MJnumcolumn_m
        WW1 = [WW1; MJtablesuccess{1,i}{1}(:,1)];
        WW2 = [WW2; MJtablesuccess{1,i}{1}(:,2)];
    end
else
    WW1=[];
    WW2=[];
    for i = 1:MJnumcolumn_m
        WW1 = [WW1; MJtablesuccess{1,i}{1}(:,1)];
        WW2 = [WW2; MJtablesuccess{1,i}{1}(:,2)];
    end  
end



%% Plot the seeding density graph of the ROI 

    histWW = [WW1,WW2];
    NN = hist3(histWW,'Ctrs',{(min(V1):grid_val:max(V1)) min(V2):grid_val:max(V2)});
    HistoAll = figure('Name', 'Seeding Histogram')
    hist3(histWW,'Ctrs',{(min(V1):grid_val:max(V1)) min(V2):grid_val:max(V2)})
    hold on 
    hm = mesh(X,Y,X*0);  
    hp = plot([x],[y],'r-') 
    set(hm,'EdgeColor','k')  
    set(hp,'LineWidth',2)  
    set(gca,'Visible','on')
    hold off
       
% % Scatter Plot to help visualise the starting points (blue) of the tracers 
% % within the region of interest (Not green)
% 
% %Uncomment from here \/
% 
% for col = MJnumcolumn_m2(1):MJnumcolumn_m2(2)
%         
%     orthopoints = figure('Name','Orthopoints')
%     hold on
% 
%     scatter(background(:,1),background(:,2),1,'green')  
% 
%     scatter(MJtablesuccess{1,col}{1}(:,1),MJtablesuccess{1,col}{1}(:,2),1,'blue')
% 
%     scatter(MJtablefailed{1,col}{1}(:,1),MJtablefailed{1,col}{1}(:,2),3,'black')
% 
%     scatter(MJtablefailed{1,col}{1}(:,3),MJtablefailed{1,col}{1}(:,4),3,'red')
%     
% end
