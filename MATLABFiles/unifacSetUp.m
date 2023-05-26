function [paramsPure] = unifacSetUp(compArray) 
%File is adapted from Elliott and Lira (https://sourceforge.net/projects/chethermo/)

%This function extracts the parameters of interest in you system from a
%large table which contains a larger variety of parameters
%Meant to be used in conjunction with unifac.m
%This table (unifacAijLLE.mat) was compiled from known parameters by
%Elliott & Lira (2012)

%compArray is an array(nGroups,nComp+1) where the value aij tells you the nu
%value of group i in component j-1 and the first column tells you the name
%of the group. See below for an example array

% For EOPO-DAB
% compArray = {
%           'OH'        2       0       0;
%           'CH2'       205     0       0;
%           'CH2O'      255     0       0;
%           'CH'        52      0       0;
%           'CH3'       52      0       0;
%           'AC'        0       2       0;
%           'ACH'       0       6       0;
%           'ACNH2'     0       4       0;
%           'H2O'       0       0       1;
%           };


%Want to return a series of arrays for the unifac function to use
%Main goal is to extract relevant parameters from a table containing all
%R, Q and interaction parameters for all groups of interest

load 'unifacAijLLE.mat'; %Table which contains all our parameters, compiled by Elliott & Lira, 2012

%Initialize some reference variables
nHeader = 2; %Tells you number of extraneous rows on top
nDataRows = size(unifacAijLLE,1);
nComp = size(compArray,2)-1;
nGroups = size(compArray,1);

%extract numerical portion of compArray
compNArray = reshape(cat(1,compArray{:,2:nComp+1}),nGroups,nComp);

%Setup variables to store our parameters of interest
r = zeros(1,nComp);
q = zeros(1,nComp);
R = zeros(1,nGroups);
Q = zeros(1,nGroups);
PsiWaterK = zeros(1,nGroups);
PsiKWater = zeros(1,nGroups);
aij = zeros(nGroups,nGroups);

%Create variable to keep track of if the group name was found, and to keep
%track of which row that group corresponds to
found = ones(1,nGroups);
groupRow = zeros(1,nGroups);

%Search table for groups of interest, and save our parameters
for group=1:nGroups
    %Loop through all groups
    for dataRow=nHeader+1:nDataRows
        %Loop through all data rows in the table
        if(strcmp(compArray{group,1},unifacAijLLE{dataRow,2})==1)
            %If the group name in compArray is the same as in the array
            found(group)=0;%register that this group has been found
            groupRow(group)=dataRow;%save what row corresponds to the group
            
            %store group volume and surface area
            R(group) = unifacAijLLE{dataRow,3};
            Q(group) = unifacAijLLE{dataRow,4};
            PsiWaterK(group) = unifacAijLLE{11,dataRow+2};
            PsiKWater(group) = unifacAijLLE{dataRow,13};
            
            %update r and q values for components
            r = r + compNArray(group,:)*R(group); 
            q = q + compNArray(group,:)*Q(group);
            break;%Group is found, so no longer need to keep looping
        end
    end
end

if(sum(found)>0) 
    %If one of the group names wasn't found
    for i=1:nGroups
        if(found(i)>0) %Find unfound group
            %Print error message
            sprintf('Group %s from row %d of compArray is not in database. Check input for typo.\n', char(compArray(i,1)), i)
        end 
    end 
    return; %exit program
end 

% fill aij matrix
% You do +2 for the columns, since the first two columns are for R,Q
for i = 1:nGroups
    for j = 1:nGroups
        aij(i,j) = unifacAijLLE{groupRow(i),groupRow(j)+2};
    end
end

%To calculate group mole fractions in pure component solutions, need to
%divide nu value for each group by the total number of groups in solution%
%In Xpure, each element aij should correspond to the mole fraction of 
%group i in a solution of pure component j

%kron() function creates copies of the sum array so that element division 
%can be done
Xpure = compNArray ./ kron(sum(compNArray,1),ones(nGroups,1));

%Next, need to calculate surface area fractions in pure component solutions
%Need to multiply Q value for group by the X value of that group and divide
%by total surface areas
calc_1 = kron(Q',ones(1,nComp)).*Xpure;
ThetaPure = calc_1./kron(sum(calc_1),ones(nGroups,1));

paramsPure = {r, q, R, Q, nComp, nGroups, compNArray, Xpure, ThetaPure, aij, PsiWaterK, PsiKWater};

end