close all
clear all
clc

%% READ HS image

fname_selected = "listy1_2023-03-21_09-40-03"; % name of the file
% fpath = 'D:\listy1_2023-03-21_09-40-03\capture\'; % path to file wihout name
fpath = 'C:\Users\krauz\Downloads\HS-image-example-20250107T082612Z-001\HS-image-example\listy1_2023-03-21_09-40-03\capture\'; % path to file wihout name

% %%
% 
% fname_selected = ''; % name of the file
% fpath = ''; % path to file wihout name

[Data, info] = enviread(strcat(fpath,fname_selected,'.raw'),strcat(fpath,fname_selected,'.hdr'));
[White, winfo] = enviread(strcat(fpath,'WHITEREF_',fname_selected,'.raw'),strcat(fpath,'WHITEREF_',fname_selected,'.hdr'));
[Dark, dinfo] = enviread(strcat(fpath,'DARKREF_',fname_selected,'.raw'),strcat(fpath,'DARKREF_',fname_selected,'.hdr'));

wavelengths = str2num(info.Wavelength(2:end-1));

%% Calibration

HS_calibrated = zeros(size(Data));
white_ref = mean(White, 1);
dark_ref = mean(Dark, 1);


for i = 1:size(Data, 1)
   HS_calibrated(i,:,:) = (Data(i,:,:) - dark_ref(1,:,:))./(white_ref(1,:,:) - dark_ref(1,:,:));
end

%% SHOW fake RGB image 
binning = 8; % set original binning factor

R = round(333/binning); % channel number (not wavelenght) / binning factor;
G = round(205/binning); % channel number (not wavelenght) / binning factor;
B = round(100/binning); % channel number (not wavelenght) / binning factor;

figure()
imshow(HS_calibrated(:,:,[R, G, B]),[])

%%  S-golay filtering
hsfiltered = zeros(size(HS_calibrated));
order = 2; % set order for SG filter
window = 15; % set window of spectral channels for SG filter 

for x = 1:length(HS_calibrated(:,1,1))
    for y = 1:length(HS_calibrated(1,:,1))
        hsfiltered(x, y, :) = sgolayfilt(HS_calibrated(x, y, :), order, window);
    end
end


%%  Seletcion of region of interest (ROI)

f = figure();
imshow(hsfiltered(:,:,[R, G, B]),[])
h_rect = imrect();
pos_rect = round(h_rect.getPosition());
close(f)
I_roi = hsfiltered(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)),:); % area of image

%% Caluclate mean reflectance from ROI
sz = size(I_roi);

I = reshape(I_roi, [sz(1) * sz(2), sz(3)]);
mean_ref = mean(I);
std_ref = std(I); 

%% Plot mean reflectance

figure()
plot(wavelengths, mean_ref,'b-')
hold on
plot(wavelengths, mean_ref + std_ref, 'r-')
plot(wavelengths, mean_ref - std_ref, 'r-')
axis([400 1000 0 1.2]) 

xlabel('\lambda (nm)')
ylabel('Normalized reflectance (-)')


%% Calulate some spectral measures (additional (Add on) hyperspectral toolbox for image processing toolbox necessary)
% https://www.mathworks.com/help/images/hyperspectral-image-processing.html

%select_men_ferfeltace for comaprisom 

reflectance1 = select_mean_ref(hsfiltered,R,G,B);

reflectance2 = select_mean_ref(hsfiltered,R,G,B);


sid_out = sid(reflectance1,reflectance2)
sam_out = sam(reflectance1,reflectance2)
sidsam_out = sidsam(reflectance1,reflectance2)
jmsam_out = jmsam(reflectance1,reflectance2)
ns3_out = ns3(reflectance1,reflectance2)


function mean_ref_fc = select_mean_ref(hsfiltered,R,G,B)
    f = figure();
    imshow(hsfiltered(:,:,[R, G, B]),[])
    h_rect = imrect();
    pos_rect = round(h_rect.getPosition());
    close(f)
    I_roi = hsfiltered(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)),:); % area of image
    
    sz = size(I_roi);
    I = reshape(I_roi, [sz(1) * sz(2), sz(3)]);
    mean_ref_fc = mean(I);
    % std_ref = std(I); 
end



