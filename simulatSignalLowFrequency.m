clc;
clear;
load('collectedData4700to7150');
t = reshape(1/48000:1/48000:10,480000,1);
numofFre = 8;

TEMPERATURE = 16;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;
f = 4350;
cosSignal = zeros(480000,numofFre);
sinSignal = zeros(480000,numofFre);

%8个频率下的signal
for freNum = 1:numofFre
    fre = f + freNum * 350;
    cosSignal(:,freNum) = reshape(cos(2*pi*fre*t),480000,1);
    sinSignal(:,freNum) = reshape(-sin(2*pi*fre*t),480000,1);
end

receivedSignal = R;
framelength = 1920;
deciSize = framelength/16;
DCvalueI = zeros(numofFre,1);
DCvalueR = zeros(numofFre,1);
maxR = zeros(numofFre,1);
maxI = zeros(numofFre,1);
minR = zeros(numofFre,1);
minI = zeros(numofFre,1);
totDis = 250;
DVBaseband = zeros(deciSize,numofFre);
temp = zeros(deciSize,numofFre);
temp1 = zeros(deciSize,1);
temp2 = zeros(deciSize,1);
Hm = mfilt.cicdecim(16,17,3);
data1 = zeros(250,1); 
data2 = zeros(250,1);
t1 = reshape(0:deciSize-1,deciSize,1);
freCount = zeros(numofFre,1);
variance = zeros(numofFre,1);
numCount = 0;
for i = 1:250
    framesignal=receivedSignal((i-1)*framelength+1:i*framelength);
    sumxy = 0;
    sumy = 0;
    freCount = zeros(numofFre,1);
    for freNum = 1:8
        cos1 = cosSignal(((i-1)*framelength+1:i*framelength),freNum);
        sin1 = sinSignal(((i-1)*framelength+1:i*framelength),freNum);
        InPhase = framesignal.*sin1;
        Quard = framesignal.*cos1;
        
        %将两路信号用CICfilter处理
        InPhaseFi = filter(Hm,InPhase);
        QuardFi = filter(Hm,Quard);

        InPhasePro = double(InPhaseFi);
        QuardPro = double(QuardFi);

        %将两路信号用LEVD算法处理并求出DC vector
        [ESCInPhase,ESCQuard,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum)] = LEVD2(InPhasePro,QuardPro,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum));
 
        %除去信号中的static vector 并将两路信号转变成baseband 
        DVInPhase = InPhasePro - ESCInPhase; 
        DVQuard = QuardPro - ESCQuard;
        DVBaseband(:,freNum) = ReImToComp(DVInPhase, DVQuard);

        % 求出相位同时根据threshold继续对DC vector进行处理

        [ph,DCvalueR(freNum),DCvalueI(freNum),freCount(freNum)] = DCprocess(DVBaseband(:,freNum),maxR(freNum),minR(freNum),DCvalueR(freNum),maxI(freNum),minI(freNum),DCvalueI(freNum),freCount(freNum));
        
        fre = f + freNum * 350;


        %linear regression
        if freCount(freNum) == 1
            for a = 1:deciSize
                temp(a,freNum) = ph(a) - ph(1); 
                temp(a,freNum) = temp(a,freNum)*SoundSpeed*1000/(2*pi*fre);
            end
            sumy = sum(temp(:,freNum))+sumy;
            sumxy = sum(temp(:,freNum).*t1)+sumxy;
            numCount = numCount + 1;
        end
    end
    
    %find ignore frequency
    if numCount == 0
        data1(i) = 0;
        data2(i) = totDis;
        continue; 
    end
    deltax = numofFre*((deciSize-1)*deciSize*(2*deciSize-1)/6-(deciSize-1)*deciSize*(deciSize-1)/4);
    delta = (sumxy-sumy*(deciSize-1)/2)/deltax*numofFre/numCount;
    temp1 = delta*t1;
    varsum = 0;
    
    %get variance of each frequency
    for freNum = 1:numofFre
        if freCount(freNum) == 1
            temp2 = temp1 - temp(:,freNum);
            variance(freNum) = sum(temp2.^2);
            varsum = varsum + variance(freNum);
        end
    end
    
    for freNum = 1:numofFre
        if variance(freNum) > varsum/numCount
            freCount(freNum) = 0;    
        end
    end
    
    sumxy = 0;
    sumy = 0;
    numCount = 0;
    
    %Caluculate the distance based on the frequency after removing the ignore frequency
    for freNum = 1:numofFre
       if freCount(freNum) == 1
           sumy = sum(temp(:,freNum))+sumy;
           sumxy = sum(temp(:,freNum).*t1)+sumxy;
           numCount = numCount + 1;
       end

    end
    if numCount == 0
       data1(i) = 0;
       data2(i) = totDis;
       continue;
    end
    delta = (sumxy-sumy*(deciSize-1)/2)/deltax*numofFre/numCount;
    distance = -delta*deciSize/2;    
    data1(i) = distance; 
    totDis = totDis + distance;
    data2(i) = totDis;

end
x = 1/25:1/25:10;
plot(x,data2,'r');
xlabel('time , s');
ylabel('distance , mm');
