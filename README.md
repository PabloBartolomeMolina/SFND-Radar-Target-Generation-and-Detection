# SFND-Radar-Target-Generation-and-Detection
Small project in MATLAB for Udacity Nanodegree in Sensor Fusion in Matlab

# FMCW Waveform Design
By using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), 
chirp time (Tchirp) and slope of the chirp. Recall that the initial parameters, some already 
given and some others fixed by me according to the specified ranges are the following ones:

% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 70 m/s

%speed of light = 3e8
v_light = 3e8;
% range resolution
r_res = 1;
% max range
r_max = 200;

% User Defined Range and Velocity of target
% define the target's initial position and velocity. Note : Velocity
% remains contant
p0 = 50;
v0 = -30;

With this data, the Bandwidth (B), the chirp time (Tchirp) and the slope are calculated:
B = (v_light) / (2*r_res);
Tchirp = (5.5*2*r_max) / (v_light);
slope = B / Tchirp;
disp(slope)

# Simulation Loop
The following code simulates the target movement and calculates the beat of mixed signal for every 
timestamp.
%For each time stamp update the Range of the Target for constant 
velocity. 
 r_t(i) = p0 + v0*t(i);
 td(i) = (2*r_t(i))/ (v_light);
 
 %For each time sample we need update the transmitted and
 %received signal. 
 Tx(i) = cos(2*pi*(fc*t(i)+slope*t(i)^2/2));
 Rx(i) = cos(2*pi*(fc*(t(i)-td(i))+slope*((t(i)-td(i))^2)/2));
 
 %Now by mixing the Transmit and Receive generate the beat signal
 %This is done by element wise matrix multiplication of Transmit and
 %Receiver Signal
 Mix(i) = Tx(i)*Rx(i);
 
 # Range FFT (1st FFT)
 The following code implements the Range FFT on the Beat or Mixed Signal and plots the result.
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the 
size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix,[Nr, Nd])
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
sig_fft = fft(Mix,Nr)./Nr;
% Take the absolute value of FFT output
sig_fft = abs(sig_fft);
% Output of FFT is double sided signal, but we are interested in only one 
side of the spectrum.
% Hence we throw out half of the samples.
sig_fft = sig_fft(1:(Nr/2));
%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)
% plot FFT output 
plot(sig_fft); %grid minor 
axis ([0 200 0 1]);
yticks(0:0.1:1);
xlabel('Measured Range');

# doppler FFT (2st FFT)
Mix=reshape(Mix,[Nr,Nd]);
% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);
% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;
%use the surf function to plot the output of 2DFFT and to show axis in 
both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

# 2D CFAR
Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map.
Determine the number of Training cells for each dimension. Similarly, pick the number of guard cells.
This is done with a try and error processes to get the result as similar as possible to the one seen in 
the walkthrough in the end. The effects of selecting training and guard cells zones too big or too 
small are seen, as well as the effect of offsets lower than 1.2 (not good at all) or bigger than 1.75 (we 
loose a lot of information in the final graph).
%Select the number of Training Cells in both the dimensions.
Tr = 8;
Td = 4;
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 6;
Gd = 3;
% offset the threshold by SNR value in dB
offset=1.25;
Slide the cell under test across the complete matrix. Make sure the CUT has margin for Training and 
Guard cells from the edges.
for i = 1:(Nr/2-(2*Gr+2*Tr+1))
 for j = 1:(Nd-(2*Gd+2*Td+1))
 …
 end
end
The code inside this double loop is explained in the following lines.
For every iteration sum the signal level within all the training cells. Convert the value from 
logarithmic to linear using db2pow function in order to perform the sum.
noise_level = zeros(1,1);
 for x = i:(i+2*Tr+2*Gr) 
 noise_level = [noise_level, db2pow(RDM(x,(j:(j+ 2*Td+2*Gd))))];
 end 
 sum_cell = sum(noise_level);
 noise_level = zeros(1,1);
 for x = (i+Tr):(i+Tr+2*Gr) 
 noise_level = [noise_level, db2pow(RDM(x,(j+Td):(j+Td+2*Gd)))];
 end 
 sum_guard = sum(noise_level);
 sum_train = sum_cell - sum_guard;
Average the sum values for all of the used training cells and then convert this value back to 
logarithmic scale using pow2db. Finally the offset is multiplied to this value to get the proper value of 
the threshold.
threshold = pow2db(sum_train/Tcell)*offset;Next, compare the signal under CUT against this threshold. If the CUT level > threshold assign it a 
binary value (0 / 1).
signal = RDM(i+Tr+Gr, j+Td+Gd);
 if (signal < threshold)
 signal = 0;
 else
 signal = 1;
 end 
 CFAR(i+Tr+Gr, j+Td+Gd) = signal;
To keep the map size same as it was before CFAR, all the non-thresholded cells of the edges are 
equated to 0.
for i = 1:(Nr/2)
 for j = 1:Nd
 if (i > (Tr+Gr))& (i < (Nr/2-(Tr+Gr))) & (j > (Td+Gd)) & (j < (Nd-
(Td+Gd)))
 continue
 end
 CFAR(i,j) = 0;
 end
end
Selection of Training, Guard cells and offset. Training, Guard cells and offset are selected by 
increasing and decreasing to match the image seen in the walkthrough. The result is presented in the 
following image.
