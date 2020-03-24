%%
clc 
clear all

% --Creating the network-- %

% Declaring Parameters
trPwr = 0.1; % Transmitting power, assumed to be same for all transmitters 
thermalNoise = 3.34 *10^-9; % Assumed to be same for all receivers
W = 3*10 ^7; % Frequency bandwidth.
gamma = 0.33;

% Declaring the nodes 
nodesM=[10,13.5;20,13.5;30,13.5;40,13.5];
baseStation = [20,30];
nNodes =length(nodesM(:,1)); 

% Plot the nodes 
figure(1)
plot(nodesM(:,1),nodesM(:,2),'o') 
axis([0 50 0 40]); 
grid on 
hold on
plot(baseStation(1),baseStation(2),'o','MarkerSize',10)
xlabel('Horizontal Position (m)');
ylabel('Vertical Position (m)');

%Total demand for each node
demand = [1397,3328,1708,1499]*10^6; %mBits

%Calculate the distance for each node to basestation
for i = 1:nNodes
    dist(i) = sqrt(((nodesM(i,1)-baseStation(1))^2)+((nodesM(i,2)-baseStation(2))^2)); 
end

% --Propagation model-- %

% Distance for the propagation model
propd = 1:(max(dist)+5);

%Parameters for the Propagation model, Friis
cKonst = 300*1000000;
freq = 30000000;
Gtx = 10;
Grx = 6.3;
PtxdBm = 43; % 25watt ish
lamb = cKonst / freq;

%Calculate the propagation for each distance
for  i=1:length(dist)
    tmp = (PtxdBm + Gtx + Grx) + 20 * log10(lamb/(4*pi*dist(i)));
    p(i) = 10^((tmp-30)/10);
end

for  i=1:length(propd)
    tmp  = (PtxdBm + Gtx + Grx) + 20 * log10(lamb/(4*pi*propd(i)));
    prop(i) = 10^((tmp-30)/10);
end

% %Plot the propagation model
figure(2)
plot(propd,prop)

% --Calculating the channelwidth for each user -- %

%Declaring parameters
N=nNodes;
f = ones(1,N);
x_assumed = ones(1,N)*1000; % Initial channel-values
yk = p'; %Received power for each user
yi = ones(1,N)*gamma; % Threshold array

% Creating a and b matrix for the opt-problem
for i = 1:N
     a(i,i) = x_assumed(i)*yk(i);
     b(i) =  gamma*(sum(x_assumed(:).*yk(:))-x_assumed(i)*yk(i));
end
% Gain
gain = linprog(f,-a,-b,[],[],zeros(1,N));

%Calculates the SINR and transform it to power in watt
for i = 1:nNodes 
        upper = p(i)*gain(i); 
        lower = thermalNoise + sum(p(:).*gain(:) - (p(i)*gain(i))); 
        SINR(i) = 1 * (10^((upper/lower)/10))/1000; 
end

% Basestation capacity
bsCapacity = 100*10^6;

% Calculate the shannon capacity for each user, floor it to int numbers of
% bits
for i = 1:nNodes 
        c(i) = floor(W * log2(1+SINR(i)));
end

% Declaring numbers for the main opt-problem
transferRate = c
nLengths = length(transferRate);

for  i = 1:nNodes
    quantity(i,1) = round(demand(i)/transferRate(i)); 
    %quantity(i,1) = round(floor(bsCapacity/transferRate(i))*mupp(i));
end

quan = quantity



patterns = diag(floor(bsCapacity./transferRate));
%patterns = diag(ones(1,nLengths));
nPatterns = size(patterns,2);

lb2 = zeros(nLengths,1);
A2 = transferRate;
b2 = bsCapacity;

lpopts = optimoptions('linprog','Display','off');
ipopts = optimoptions('intlinprog',lpopts);

reducedCost = -Inf;
reducedCostTolerance = -0.0001;

exitflag = 1;

% while reducedCost < reducedCostTolerance && exitflag > 0
       lb = zeros(nPatterns,1);
%     f = lb + 1;
       A = -patterns;
        b = -quantity;
%     
%     [values,nSeconds,exitflag,~,lambda] = linprog(f,A,b,[],[],lb,[],lpopts); 
%  
%     if exitflag > 0
%         fprintf('Using %g seconds\n',nSeconds);
%         % Now generate a new pattern, if possible
%         f2 = -lambda.ineqlin;
%         [values,reducedCost,pexitflag] = intlinprog(f2,1:nLengths,A2,b2,[],[],lb2,[],ipopts);
%         reducedCost = 1 + reducedCost; % continue if this reducedCost is negative
%         newpattern = round(values);
%         if pexitflag > 0 && reducedCost < reducedCostTolerance
%             patterns = [patterns newpattern];
%             nPatterns = nPatterns + 1;
%         end
%     end
% end

cnt11=1;
if exitflag <= 0 
    disp('Error in column generation phase')
else
    [values,secondsUsed,exitflag] = intlinprog(f,1:length(lb),A,b,[],[],lb,[],[],ipopts);
    if exitflag > 0
        values = round(values);
        secondsUsed = round(secondsUsed);
        fprintf('Optimal solution uses %g seconds\n\n', secondsUsed);
        totalwaste = sum((patterns*values - quantity).*transferRate'); % waste due to overproduction
        
        for j = 1:size(values)
            if values(j) > 0
                totThroughput=0;
                fprintf('In %g seconds\n\n',values(j));
                for w = 1:size(patterns,1)
                    if patterns(w,j) > 0
                        fprintf('User %d sends %d packets with size %d bit\n', w,patterns(w,j),transferRate(w));
                        waaste(cnt11,1)=values(j);
                        totThroughput = totThroughput + (patterns(w,j)*transferRate(w));
                    end
                end
                wastej = bsCapacity - dot(patterns(:,j),transferRate); % waste due to pattern inefficiency
                totalwaste = totalwaste + wastej;
            fprintf('    Waste of this pattern is %g\n\n', wastej);
           waaste(cnt11,2)=totThroughput;
           waaste(cnt11,3)=bsCapacity-totThroughput; 
           cnt11 = cnt11+1;
            end
        end
        fprintf('Total waste in this problem is %g.\n',totalwaste);
    else 
        disp('Error in final optimization')
    end
end

B=sortrows(waaste,-2);

tt=1;
for i = 1:length(waaste(:,1)) 
    wasteToPlot(tt:tt+waaste(i,1)-1,2)=B(i,2);
    wasteToPlot(tt:tt+waaste(i,1)-1,3)=B(i,3);
    tt = tt+waaste(i,1);
end
wasteToPlot(:,1)=1:length(wasteToPlot(:,2));
figure(3)
plot(wasteToPlot(:,1),wasteToPlot(:,2))
xlabel('Second');
ylabel('Throughput');

figure(4)
plot(wasteToPlot(:,1),wasteToPlot(:,3))
xlabel('Second');
ylabel('Waste');
