load('collectedData1900to7150.mat');
t = reshape(1/48000:1/48000:10,480000,1);
numofFre = 16;
framelength = 1920;
deciSize = framelength/16;
DCvalueI = zeros(numofFre,1);
DCvalueR = zeros(numofFre,1);
maxR = zeros(numofFre,1);
maxI = zeros(numofFre,1);
minR = zeros(numofFre,1);
minI = zeros(numofFre,1);
totDis = 0;
baseband = zeros(16,1);
DVBaseband = zeros(deciSize,numofFre);
temp = zeros(deciSize,numofFre);
temp1 = zeros(deciSize,1);
temp2 = zeros(deciSize,1);
Hm = mfilt.cicdecim(16,17,3);
data1 = zeros(250,1); 
data2 = zeros(250,1);
t1 = reshape(0:deciSize-1,deciSize,1);
freCount = zeros(numofFre,1);
var = zeros(numofFre,1);
numCount = 0;
baseband1 = zeros(250,numofFre);


TEMPERATURE = 16;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;
f = 1550;
receivedSignal = R;
cosSignal = zeros(480000,numofFre);
sinSignal = zeros(480000,numofFre);

%16个频率下的signal
for freNum = 1:numofFre
    fre = f + freNum * 350;
    cosSignal(:,freNum) = reshape(cos(2*pi*fre*t),480000,1);
    sinSignal(:,freNum) = reshape(-sin(2*pi*fre*t),480000,1);
end

for i = 3:250
    framesignal=receivedSignal((i-1)*framelength+1:i*framelength);

    for freNum = 1:numofFre
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
        
        [ph,DCvalueR(freNum),DCvalueI(freNum),freCount(freNum)] = DCprocess(DVBaseband(:,freNum),maxR(freNum),minR(freNum),DCvalueR(freNum),maxI(freNum),minI(freNum),DCvalueI(freNum),freCount(freNum));
       
        baseband1(i,freNum) = sum(DVBaseband(120,freNum));
    

    end
    
    baseband = reshape(baseband1(i,:),16,1);
    ift = dsp.IFFT;    
    y = ift(baseband);

    [~, indy] = max(abs(y));

    data1(i) = (indy-1)*0.0615;

end
x = 1/25:1/25:10;

plot(x,data1,'r');
xlabel('time , s');
ylabel('distance , m');
legend('预期距离' , 'IDFT测量距离');

load train;
sound(y,Fs);