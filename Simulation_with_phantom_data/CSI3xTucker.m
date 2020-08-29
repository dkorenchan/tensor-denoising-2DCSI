function CSI3xTucker
% CSInoising.m - Written by D. Korenchan, based upon MATLAB code by C. Baligand
% and C. Najac
%
% Reads in 2D CSI spectra (using SIVIC systems), adds in Gaussian noise to
% all FIDs and then denoises using the tensor value decomposition
% algorithm. The noised and denoised data are plotted in an interactive
% GUI, which enables the user to customize denoising settings and run
% optimization/characterization tests.
%      
%
% USAGE:
% 1.    Change all values to adjust prior to running near start of this function
% 2.    When prompted, specify .fid file containing CSI data
%
% CHILD SCRIPTS:    load_fids.m, load_procpar.m,
%                   scripts from SIVIC software
%
% UPDATES:
%    
%   5/22/20:    Ported over from CSInoisingmod.m. Added extra steps in
%               noiseDenoise: after denoising, undo apodization by
%               multiplication in k-space, re-FFT, and denoise again
%               (undoing apodization makes the noise more uniform, enables
%               further noise removal)
%   5/23/20:    Modified calculation of noise for a desired SNR in both
%               snrTest and setSNR: use imgnoisy data as reference to set
%               SNR, rather than original img data (since apodization will
%               also affect SNR)
%   5/24/20:    Included final apodization step after 2x denoising to match 
%               apodization in original dataset (better visual comparison
%               of spectra, more accurate integral quantification). Also
%               changed setSNR back to estimate SNR from img rather than
%               imgnoisy, b/c it was not robust. Also moved procCSI calls
%               to calcMaps so ensure that map calculations are always done
%               with a fresh CSI dataset
%   5/24/20:    BIG BUG FIX: Spatial grid shifting has not been being 
%               performed since code was changed to do k-space processing 
%               independent of SIVIC!! Added code to shift data in k-space 
%               based upon diffvox variable
%   5/25/20:    Added code to perform Gaussian apodization rather than
%               Lorentzian (rationale: reduce baseline broadening that may
%               affect smaller peak). Verified that it matches SIVIC
%               apodization. 
%   5/26/20:    Minor fix: SNR test would miscalculate midpoint for SNR
%               values to test if set too low; fixed by using setSNR
%               subfunction beforehand.
%   5/27/20:    Added an initial denoising step (only 1/2 truncation in
%               spectral) prior to any apodization. 2x Tucker proceeds from 
%               there as normal
%

% VALUES TO ADJUST PRIOR TO RUNNING
%
origapod = 10; %apodization to use with non-denoised data, in Hz
% sizeLR_firstpass = [3,6,6]; %denoising used the first time on data
load_path = '/Users/sf865719/Lab/Data/600_Data/GLC_phantoms/2DCSI/'; 
    %path to data
% load_path = '/Users/sf865719/Lab/Data/600_Data/pH_invivo/GLC/RCC/s_20180305_01_RCC014-GLConly';
%     %path to data
save_dir = '/Users/sf865719/Lab/MATLAB/Simulation/TVD_denoising/Results';
    %path for saving data
svk_path = '/Applications/SIVIC.app/Contents/sivic/local/bin/'; %path to SIVIC scripts
% shift_path = '/Users/sf865719/Lab/Processing_macros/csishift'; %path to csishift script
filestrref = 'fsems_mouse_axial'; %subset of reference image filename
% filestrref = 'sems_mouse_axial_ldha'; %subset of reference image filename
% filestrref = 'ref'; %set to 'ref' if loading from ref.idf
spatzfillfactor = 1; %factor for spatial zerofilling
% globshft = [-0.3,1.0]; %[mm right,mm down] shift of entire 13C data relative to 1H
%     %(determined for 11/17/15 phantom data)
% globshft = [0.4,1.4]; %[mm right,mm down] shift of entire 13C data relative to 1H
%     %(determined using 2/6/18 phantom data)
globshft = [-0.2,0.3]; %[mm right,mm down] shift of entire 13C data relative to 1H
    %(determined for 3/3/18 phantom data)    
% globshft = [-2,1.4]; %[mm right,mm down] shift of entire 13C data relative to 1H
%     %(determined specifically for 8/14/15 phantom data)
nSVspatmax = 8; %used for SV test
brukflg = false; %set to true if loading from Bruker ser file
lastflg = false; %selects final file found when searching with filestrref;
    %otherwise, specified by user
nopeflg = false; %makes things less verbose for multirun    
tol = 0.15; %allowable pH tolerance 
% sneakmask = [0  0   0   0   0   0   0   0; ...
%             0   0   0   1   1   0   0   0;...
%             0   0   0   1   1   0   0   0;...
%             0   1   1   0   1   1   0   0;...
%             0   1   1   0   1   1   0   0;...
%             0   0   0   0   0   0   0   0;...
%             0   0   0   0   0   0   0   0;...
%             0   0   0   0   0   0   0   0]; %used for 2/6/18 phantom data,
%                                             %to isolate tube voxels
% sneakmask = [0  0   0   0   0   0   0   0; ...
%             0   0   0   0   0   0   0   0;...
%             0   1   1   1   1   1   0   0;...
%             0   1   1   1   1   1   0   0;...
%             0   1   1   1   1   1   0   0;...
%             0   1   1   1   1   1   0   0;...
%             0   0   1   1   1   1   0   0;...
%             0   0   0   0   0   0   0   0]; %used for 3/5/18 in vivo RCC data,
%                                             %to isolate mouse voxels
% sneakmask = ones(8);                                            

close all
pKa = 6.17; %pKa of bicarbonate-CO2
fcorr = 9.276; %flip angle correction factor for pH imaging
spnames = {'BiC','CO2'}; %specifies names of spectral bins to define
imgnames = {'img','imgapod','imgapodnoisy','imgnoisy','imgdenoisedFINAL'};
maps.nRuns = 5; %# of times to run noising+denoising for evaluating 
    %performance (adjustable via GUI)
nSVspecmax = 8; %used for SV test (adjustable via GUI)
colors = [1 0 0;0 1 0;0 1 1;1 1 0;1 0 1;0 0 1];
lims = struct;
home = pwd;

% % Check to see if previous images are an input. If so, load images and skip 
% % to ROI processing. Otherwise, continue with function as normal.
% %
% if nargin > 0 && exist('imageprev','var') == 1 && isfield(imageprev,'ref')
%     prompt = {'Previous images detected as input. Use?'};
%     choices = {'Yes' 'No'};
%     answer2 = listdlg('ListString' , choices , ...
%         'SelectionMode' , 'single' , ...
%         'ListSize' , [200 30] , ...
%         'PromptString' , prompt);
%     if answer2 == 1 
%         disp('Loading images from input data...');
%         images = imageprev;
%         loadflg = false;
%     else
%         disp('Images will not be loaded from previous data. Continuing...');
%         loadflg = true;
%     end
% else
    loadflg = true;
% end

%% DATA PROCESSING: Image Load

% Prompt user to specify location of .fid file containing 2D CSI data  
% as well as folder containing 1H reference .idf and .real 
% files. Load reference image data into structure ref
%

if loadflg %only perform next section if images not loaded from previous input
    while 1
        cd(load_path);
        if brukflg %load from ser file - NOT WORKING YET
            cd(pathname);
            [filename , pathname] = uigetfile('*' , 'Specify 2D CSI ser file to open'); 
                %file and path names of raw data
            if isequal(filename , 0) || isequal(pathname , 0) 
                cd(home);
                error('File not specified. Aborting function.')
            end
            % Load procpar info from 2D CSI .fid
            images.CSI.fov = [32 32 10];
            %FOV, in mm: 1st PE is x, 2nd PE is y
            images.CSI.sw = 2000; %spectral width, in Hz
            images.CSI.npx = [8 8 1]; %# of pixels
            images.CSI.off = [globshft 0];
            images.CSI.ps = images.CSI.fov ./ images.CSI.npx; %pixel size, in mm
            images.CSI.np = 128;
            images.CSI.sr = images.CSI.sw / images.CSI.np; %spectral resolution, in Hz
        else
            [filename , pathname] = uigetfile('*.fid' , 'Specify 2D CSI .fid file to open'); 
            %file and path names of raw data
            if isequal(filename , 0) || isequal(pathname , 0) 
                cd(home);
                error('File not specified. Aborting function.')
            end
            strings = split(pathname,'/');
            savename = strings{end-1};
            % Load procpar info from 2D CSI .fid directory
            cd([pathname filename]);
            info = load_procpar('procpar');
            images.CSI.fov = [info.lpe info.lpe2 (info.thk / 10)] * 10;
            %FOV, in mm: 1st PE is x, 2nd PE is y
            images.CSI.sw = info.sw; %spectral width, in Hz
            images.CSI.npx = [info.nv info.nv2 1]; %# of pixels
            images.CSI.np = info.np / 2; %# of spectral points
            images.CSI.ps = images.CSI.fov ./ images.CSI.npx; %pixel size, in mm
            images.CSI.sr = images.CSI.sw / images.CSI.np;
            %spectral resolution, in Hz
            images.CSI.off = [info.ppe+globshft(1) info.ppe2+globshft(2) ...
                info.pss0];
        end
        % Display acquisition parameters
        disp('2D CSI acquisition parameters:')
        disp(['Spatial: ' num2str(images.CSI.fov(1)) ' x ' ...
            num2str(images.CSI.fov(1)) ' mm^2 FOV, ' ...
            num2str(images.CSI.npx(1)) ' x ' ...
            num2str(images.CSI.npx(2)) ' matrix'])
        disp(['Spectral: ' num2str(images.CSI.sw) ' Hz spectral width, ' ...
            num2str(images.CSI.np) ' complex points, ' ...
            num2str(images.CSI.sr) ' Hz spectral resolution'])
        % Load raw 2D CSI k-space data
        if ~brukflg
            cd([pathname filename]);
            system([svk_path 'svk_file_convert -i fid -o csixy -t2']);
            ksp = read_ddf_image('csixy');
            cd(home);
            ksp.img = permute(ksp.img,[3 2 1]); %make matrix y by x by f
            def = [images.CSI.npx(1)*spatzfillfactor , images.CSI.np , ...
               origapod , 0 , 0];
            procCSI(def,'img');
        end
        % Read in 1H reference image
        cd(pathname);
        if strcmp(filestrref,'ref') %check if ref.idf is specified
            disp('ref image: Loading ref.idf')
            blah = read_idf_image('ref'); %load from .idf
            images.ref.img = permute(blah.img,[2,1,3]);
            images.ref.fov = blah.idf.fov; %fov lengths, in mm
            images.ref.np = blah.idf.npix; %# of pixels
            images.ref.ps = blah.idf.pixelsize; %pixel dimensions, in mm
            images.ref.off = blah.idf.LPScenter; %offsets from center, in mm
        else %load from .fid
            fnames = dir([filestrref '*.fid']);
            filenames = cell(length(fnames),1);
            for j = 1:length(fnames) %read in file names
                filenames{j} = fnames(j).name;
            end
            if lastflg
                H1filename = filenames{end};
            else
                prompt = 'Choose file to load:';
                choices = filenames;
                answer3 = listdlg('ListString' , choices , ...
                         'SelectionMode' , 'single' , ...
                         'ListSize' , [200 50] , ...
                         'PromptString' , prompt);
                H1filename = filenames{answer3};
            end
            disp(['ref image: Loading ' H1filename])
            cd(H1filename);
            H1info = load_procpar('procpar');
            %determine which parameters relate to x, y, z; load parameters
            f_id = [{H1info.dimX},{H1info.dimY},{H1info.dimZ}]; %parameter names for xyz-FOV lengths
            o_id = [{H1info.posX},{H1info.posY},{H1info.posZ}]; %parameter names for xyz-offsets
            paramload(f_id,o_id,'ref')
            cd('..')
            [re,im] = load_fids(strtok(H1filename,'.'));
            ksp.ref.raw = refimgorder(re,im,images.ref.np,images.ref.etl);
            for j = 1:images.ref.np(3)
                images.ref.img(:,:,j) = fftshift(fftn(fftshift(squeeze(ksp.ref.raw(:,:,j)))));
            end
            images.ref.img = abs(images.ref.img);
        end
        cd(home);
        break;
    end
    cd(home);

    % Determine vectors for displaying all images within same FOV + scaling
    % factors
    %
    images.plot.Hxbnds = [-1*images.ref.fov(1),images.ref.fov(1)] ...
    / 2 + images.ref.off(1); %for 1H, axial
    images.plot.Hybnds = [-1*images.ref.fov(2),images.ref.fov(2)] ...
    / 2 + images.ref.off(2);
    images.plot.Hzbnds = [-1*images.ref.fov(3),images.ref.fov(3)] ...
    / 2 + images.ref.off(3);
    %C13 plot vectors moved to procCSI function, to be updated properly

end

% Calculate global variables to be used later for image display
%
pixrat = images.CSI.npx(1) / (images.CSI.npx(1) + 1);
% offHa = (images.ref.off(3) - images.CSI.off(3)) / ...
%     images.ref.ps(3); %difference in slice stack offsets, 1H vs 2D CSI
% nsHa = round(images.CSI.ps(3) / images.ref.ps(3) / 2) * 2; 
%     %number of 1H slices in 13C slice
% stHa = round((images.ref.np(3) / 2) - (nsHa / 2 - 1) - offHa);%-1;
% fnHa = round((images.ref.np(3) / 2) + (nsHa / 2) - offHa);%-1;

%% DATA PROCESSING: Add Gaussian Noise, Denoise
%
noiseamp = max(reshape(abs(images.CSI.img),1,[])) / 18; %amplitude of noise
sizeLR = [3, 6, 6, 8, 8]; %number of singular values to KEEP for each dimension
nSVspec = 0;
nSVspat = 0;

% apod = origapod; %so noiseDenoise uses specified apodization initially
% noiseDenoise;


%% DATA PROCESSING: Prepare to Compute Maps

% Plot all spectra, prompt user to specify spectral bins of interest and 
% noise (NOTE: spectral/noise bins are only specified once!)
%
% defbinsname = 'imgdenoised';
defbinsname = 'img';

figure; hold on
for i = 1:images.CSI.npx(1)
    for j = 1:images.CSI.npx(2)
        nspec = squeeze(abs(images.CSI.(defbinsname)(i,j,:)) ...
            / max(reshape(abs(images.CSI.(defbinsname)),1,[])));
        plot(nspec)
    end
end
pause(1); %otherwise, won't plot spectra for some reason...
for i = 1:length(spnames)
    spname = spnames{i};
    title([spname ': Select spectral region'])
    nb = getrect;
    spbin.(spname) = round(nb(1)):round(nb(1) + nb(3));
    yvec = [0 1];
    xvec = [1 1] * spbin.(spname)(1);
    line(xvec,yvec,'Color',colors(i,:),'LineStyle','--')
    xvec = [1 1] * spbin.(spname)(end);
    line(xvec,yvec,'Color',colors(i,:),'LineStyle','--')
end
title('Select noise region')
nb = getrect;
spbin.noise = round(nb(1)):round(nb(1) + nb(3));
xvec = [1 1] * spbin.noise(1);
line(xvec,yvec,'Color',colors(3,:),'LineStyle','--')
xvec = [1 1] * spbin.noise(end);
line(xvec,yvec,'Color',colors(3,:),'LineStyle','--')
    

%% DATA PROCESSING: Compute, Plot Maps

% Plot interactive overlay
%
nsax = size(images.ref.img,3);
axsl = round(nsax / 2);
diffvox = zeros(2,1); %used for voxel shifting
zoomfac = 0; %factor used in zooming (positive zooms in; negative zooms out)
finishflg = false; %becomes true when finished with figure 
multdispflg = true; %toggles selected voxel outlines on and off
brkflg = false; %used to break out of SV or SNR test
saveflg = false; %used to indicate if SV/SNR test results saved
scrsz = get(groot,'ScreenSize');
sppi = get(groot,'ScreenPixelsPerInch');

ovlynames = {['Integral_' spnames{1} '_original'],...
    ['Integral_' spnames{2} '_original'],['LW_' spnames{1} '_original'],...
    ['LW_' spnames{2} '_original'],['LW_' spnames{1} '_apodized'],...
    ['LW_' spnames{2} '_apodized'],'pH_original','pH_denoised'...
    ['avgUpSNR_' spnames{1}],['avgUpSNR_' spnames{2}],'avgdpH',...
    'pctgoodSNR','pctgoodpHgoodSNR','pctrescue','pctlost'};
ovlyname = ovlynames{1};
otp = 0.3; %transparency of 13C overlays
zoomfac = 0;

% Display 1H image
%
disp('Displaying images...')
voxfig = figure(100);
set(voxfig,'Position',[0,scrsz(4)*1/6,scrsz(3),scrsz(4)*5/6]*sppi,...
    'Units','inches');
dispvox = false(images.CSI.npx(1)); %used to display spectra from a single voxel
selvox = false(images.CSI.npx(1)); %used to select multiple voxels for SV/SNR tests
spax = [];

cutoff = 3; %SNR cutoff for masking
SNRset = cutoff; %desired SNR for given selected voxel
nvox = 3; %# of voxels to select for SV/SNR test

apod = origapod; %apodization, in Hz (adjustable via GUI)
recalcAll; %in order to process original dataset
multCalc; plotAxOvly;

% Create uicontrols in figure for selecting voxels, viewing spectra
%
bg1 = uibuttongroup('Position',[0 0 1 .15]);
bg2 = uibuttongroup('Position',[0 .15 .1 .85]);
%toggle button for voxel selection
tb = uicontrol(bg1,'Style','pushbutton','Position',[0 60 160 30],...
    'String','Turn individual voxel display off','Callback',@dispVoxOnOff);
%slice selection
uicontrol(bg1,'Style','text','Position',[0 32 160 20],...
    'String','Active Slice, Axial');
uicontrol(bg1,'Style','slider','Min',1,'Max',nsax,'Value',axsl,...
    'SliderStep',[1/nsax 3/nsax],'Position',[0 0 140 30],'Callback',@axSelect);
%select active image data
uicontrol(bg1,'Style','text','Position',[160 32 160 20],...
    'String','Active Overlay Data');
uicontrol(bg1,'Style','popupmenu','Position',[160 20 160 10],...
    'String',ovlynames,'Callback',@selOvly);
%set noise level, denoising parameters
uicontrol(bg1,'Style','text','Position',[360 70 160 20],...
    'String','Desired SNR (displayed voxel):');
snrs = uicontrol(bg1,'Style','edit','Position',[400 50 80 20],...
    'String',num2str(SNRset),'Callback',@setSNR);
uicontrol(bg1,'Style','text','Position',[400 30 80 20],...
    'String','Noise Level:');
na = uicontrol(bg1,'Style','edit','Position',[400 10 80 20],...
    'String',num2str(noiseamp),'Callback',@setNoise);
uicontrol(bg1,'Style','text','Position',[500 45 80 30],...
    'String',['Spectral SVs Kept: (1-' num2str(nSVspec) ')']);
scsv = uicontrol(bg1,'Style','edit','Position',[500 10 80 20],...
    'String',num2str(sizeLR(1)),'Callback',@setSpecSV);
uicontrol(bg1,'Style','text','Position',[600 45 80 30],...
    'String',['Spatial SVs Kept: (1-' num2str(nSVspat) ')']);
stsv = uicontrol(bg1,'Style','edit','Position',[600 10 80 20],...
    'String',num2str(sizeLR(2)),'Callback',@setSpatSV);
uicontrol(bg1,'Style','text','Position',[700 45 80 30],...
    'String','SNR Cutoff:');
snrc = uicontrol(bg1,'Style','edit','Position',[700 10 80 20],...
    'String',num2str(cutoff),'Callback',@setCutoff);
uicontrol(bg1,'Style','text','Position',[800 45 80 30],...
    'String','# of Runs:');
nr = uicontrol(bg1,'Style','edit','Position',[800 10 80 20],...
    'String',num2str(maps.nRuns),'Callback',@setRuns);
%noising, denoising, voxel test buttons
% uicontrol(bg1,'Style','pushbutton','Position',[900 0 140 40],...
%     'String','NOISE AND DENOISE','Callback',@recalcAll);
mr = uicontrol(bg1,'Style','pushbutton','Position',[900 60 140 30],...
    'String',['MULTIRUN (x' num2str(maps.nRuns,'%i') ')'],...
    'Callback',@multCalc);
uicontrol(bg1,'Style','pushbutton','Position',[900 30 140 30],...
    'String','SV TEST (selected voxels)','Callback',@svTest);
uicontrol(bg1,'Style','pushbutton','Position',[900 0 140 30],...
    'String','SNR TEST (selected voxels)','Callback',@snrTest);
uicontrol(bg1,'Style','text','Position',[1050 60 80 30],...
    'String','SV test: spectral max')
ssvsm = uicontrol(bg1,'Style','edit','Position',[1050 30 80 20],...
    'String',num2str(nSVspecmax),'Callback',@setSVspecmax);
uicontrol(bg1,'Style','pushbutton','Position',[1050 0 120 30],...
    'String','ABORT TEST','Callback',@abortTest);
uicontrol(bg1,'Style','checkbox','Position',[1210 70 120 20],...
    'String','Save test results','Callback',@setSave);
%finish button
uicontrol(bg1,'Style','pushbutton','Position',[1210 0 100 60],...
    'String','FINISH','Callback',@finishFig);
%apodization + voxel shifting button
uicontrol(bg2,'Style','text','Position',[10 500 80 30],...
    'String','Spectral apodization (Hz):')
uicontrol(bg2,'Style','edit','Position',[10 470 80 20],...
    'String',num2str(apod),'Callback',@apodData);
uicontrol(bg2,'Style','pushbutton','Position',[0 430 120 30],...
    'String','Grid shift CSI data','Callback',@gridShift);
%grid/voxel toggle buttons + slider for voxel selection overlay
uicontrol(bg2,'Style','text','Position',[0 360 120 40],...
    'String','Overlay Transparency:');
uicontrol(bg2,'Style','slider','Min',0,'Max',1,'Value',otp,...
    'SliderStep',[0.1 0.3],'Position',[0 340 120 30],'Callback',@setOvlyTransparency);
uicontrol(bg2,'Style','slider','Min',-0.5,'Max',0.5,'Value',zoomfac,...
    'SliderStep',[0.01 0.1],'Position',[0 300 120 30],'Callback',@setZoom);
uicontrol(bg2,'Style','pushbutton','Position',[0 270 120 30],...
    'String','Reset zoom','Callback',@resetZoom);
%multiple voxel selection tools
uicontrol(bg2,'Style','pushbutton','Position',[0 210 120 30],...
    'String','Select multiple voxels','Callback',@selMultVox);
uicontrol(bg2,'Style','text','Position',[0 180 120 30],...
    'String','Number of voxels to select:');
uicontrol(bg2,'Style','edit','Position',[20 160 80 20],...
    'String',num2str(nvox),'Callback',@setMultVoxNum);
tbmv = uicontrol(bg2,'Style','pushbutton','Position',[0 120 120 30],...
    'String','Turn voxel outlines off','Callback',@dispMultVoxOnOff);
%update text when running tests, procedures
ut = uicontrol(bg2,'Style','text','Position',[0 20 120 80],...
    'String','');

dispVox;
waitfor(voxfig);


% %% DATA PROCESSING: Save Data 
% % Open browser window to save data. Variables saved are images and maps.
% %
% while 1
%     cd(save_dir)
%     [sfile , spath] = uiputfile('*.mat' , 'Save Workspace Variables As');
%     if isequal(sfile , 0) || isequal(spath , 0)
%         prompt = {'Are you sure you do not want to save your data?'};
%         choices = {'Yes' 'No'};
%         answer = listdlg('ListString' , choices , ...
%             'SelectionMode' , 'single' , ...
%             'ListSize' , [200 30] , ...
%             'PromptString' , prompt);
%         if answer == 1
%             cd(home);
%             disp('Data were not saved....');
%             break;
%         end
%     else
%         save(fullfile(spath , sfile) , 'images' , 'maps');
%         cd(home);
%         disp(['All images and maps saved in ' , fullfile(spath , sfile) , '!']);
%         break;
%     end
% end


%% INTERNAL FUNCTIONS - IMAGE PROCESSING
%
% paramload: Takes loaded procpar info, saves into images._.procpar
function paramload(f_id,o_id,dname)
images.(dname).orient = H1info.orient;
images.(dname).fov = [H1info.(f_id{1}) H1info.(f_id{2}) H1info.(f_id{3})] * 10; %fov lengths, in mm
images.(dname).np(strcmp(f_id,'lpe')) = H1info.nv; %matrix sizes in xyz
images.(dname).np(strcmp(f_id,'lpe2')) = H1info.nv2;
images.(dname).np(strcmp(f_id,'lro')) = H1info.np/2;
images.(dname).ps = images.(dname).fov ./ images.(dname).np; %pixel dimensions, in mm
images.(dname).off = [H1info.(o_id{1}) H1info.(o_id{1}) H1info.(o_id{1})]; %offsets from center, in mm
if H1info.ns > 1 %identify if multiple slices performed, to adjust parameters
    images.(dname).fov(strcmp(f_id,'lpe2')) = (H1info.thk + H1info.gap*10) * H1info.ns; %adjust slice-select FOV
    images.(dname).np(strcmp(f_id,'lpe2')) = H1info.ns;
    images.(dname).ps(strcmp(f_id,'lpe2')) = H1info.thk;
end
if isfield(H1info,'etl') %for fsems images
    images.(dname).etl = H1info.etl;
else
    images.(dname).etl = 1;
end
end

% procCSI: Processes raw CSI data: spectral deletions, adding noise, 
% apodization, zerofilling, and/or peak alignment. If noiseflg = true, add
% noise to spectra with amplitude stored in noiseamp. If shiftflg = true, 
% shift over each spectrum in order to align BiC + CO2 peaks with one 
% another (shifting involves both k-space shifting and moving points over). 
% If deleteflg = true, remove spectral regions at end of processing.
% If ncflg = true, include noise conditioning following apodization
% via noise addition.
% CURRENTLY NEITHER SPECTRAL DELETIONS OR PEAK ALIGNMENT ARE FUNCTIONAL
%
function procCSI(parsin,nname)%,shiftflg,deleteflg)
% if nargin < 2
%     shiftflg = false;
%     deleteflg = true;
% elseif nargin < 3
%     deleteflg = true;
% end
% remove deletion regions specified in spbin field (if deleteflg)
temp = fftshift(fftn(ksp.img)); %to spectrum
% if deleteflg 
%     temp(:,:,images.spbin.delete) = [];
% end
temp = ifftn(temp); %back to FID
if contains(nname,'noisy') %add white Gaussian noise to k-space
    noise = (randn(size(temp)) + 1i * randn(size(temp))) * noiseamp / 256 ...
        / sqrt(2); %complex noise matrix; need to divide by sqrt(2) in 
        %order to make noise standard deviation proportional to noiseamp!
        %(b/c real + imaginary components both amplified). Also divide by
        %256 so SNR calculations are more accurate
	temp = temp + noise;    
end
dt = ksp.ddf.dwelltime/1000;
t(1,1,:) = 0:dt:dt*(size(temp,3)-1); %time vector for apodization/shifting 
% %spectral shifting via complex exponential (if shiftflg = true)
% if shiftflg
%     ishft = round(images.shiftmap); %map of integer shifts (in indices) 
%         %that will be performed via circshift
%     kshft = images.shiftmap - ishft; %map of shifts (in indices) that will
%         %be performed using FID multiplication, ranging from +/-0.5
%     kshft = kshft * ksp.ddf.sweepwidth / ksp.ddf.specpoints; 
%         %convert kshft to units of frequency    
%     shftmat = zeros(size(temp.img));
%     for ii = 1:size(temp,1)
%         for jj = 1:size(temp,2)
%             shftmat(ii,jj,:) = exp(-1i * 2 * pi * t * kshft(ii,jj));
%             %the (-) factor seems to perform better than the (+) factor
%         end
%     end
%     temp = temp .* shftmat; 
% end
%apodization
% apdmat = repmat(exp(-t*parsin(3)*pi),ksp.ddf.npix(1), ...
%     ksp.ddf.npix(2),1); %exponential; an additional factor of 2*pi is 
%     required to match SIVIC apodization results
sigma = 1 * sqrt(2 * log(2)) / pi / parsin(3); %sigma parameter of 
    %Gaussian; an additional factor of 1/(2*pi) is required to match SIVIC 
    %apodization results
apdmat = repmat(exp(-1 * (t .^ 2 / 2 / (sigma ^ 2))),ksp.ddf.npix(1), ...
    ksp.ddf.npix(2),1); %gaussian 
temp = temp .* apdmat; %apodize
%spatial grid shifting
kmax = 1 / 2 / ksp.ddf.pixel_spacing(1); %k-max (prior to zerofilling), in 1/mm
kvals = linspace(-kmax,kmax,ksp.ddf.npix(1));
shiftx = parsin(4) * ksp.ddf.pixel_spacing(1) * ksp.ddf.npix(1) / parsin(1); 
    %spatial x-shift, in mm
shifty = parsin(5) * ksp.ddf.pixel_spacing(2) * ksp.ddf.npix(2) / parsin(1); 
    %spatial y-shift, in mm
imshft = exp(1i * 2 * pi * kvals' * shifty) * ...
    exp(1i * 2 * pi * kvals * shiftx);
imshft = repmat(imshft, [1 1 size(temp,3)]);
temp = temp .* imshft;
% %noise conditioning (to restore white Gaussian noise characteristics) -
% %REJECTED!
% if strcmp(nname,'imgapodnoisy')
%     ncmat = sqrt((1 - apdmat(1,1,end)^2) - (apdmat .^ 2));
%     noise2 = (randn(size(temp)) + 1i * randn(size(temp))) * noiseamp / 16 ...
%         / sqrt(2) .* ncmat;    
%     temp = temp + noise2;
% end
%zerofilling
pad = (parsin(1) - ksp.ddf.npix(1)) / 2;
temp = padarray(temp,[1 1 0]*pad); %zerofill spatially
pad = (parsin(2) - ksp.ddf.specpoints);
temp = padarray(temp,[0 0 1]*pad,0,'post'); %zerofill spectrally
images.CSI.(nname) = fftn(temp); %to spectra
images.CSI.(nname) = flip(images.CSI.(nname),1); 
% % shifting over points (if shiftflg = true)
% if shiftflg
%     for ii = 1:size(images.CSI.(nname),1)
%         for jj = 1:size(images.CSI.(nname),2)
%             images.CSI.(nname)(ii,jj,:) = circshift(images.CSI.(nname)(ii,jj,:),...
%                 - ishft(ii,jj)); %ishft needs to be negative for best results
%         end
%     end
% end
% Update CSI parameters based upon zerofilling + loading in image 
images.CSI.fov = [info.lpe info.lpe2 (info.thk / 10)] * 10;
%FOV, in mm: 1st PE is x, 2nd PE is y
images.CSI.sw = info.sw; %spectral width, in Hz
images.CSI.npx = [parsin(1) parsin(1) 1]; 
    %# of pixels
images.CSI.ps = images.CSI.fov ./ images.CSI.npx;
    %#pixel size, in mm
images.CSI.off = [info.ppe+globshft(1)+parsin(4)*images.CSI.ps(1) ...
    info.ppe2+globshft(2)+parsin(5)*images.CSI.ps(2) info.pss0];    
images.CSI.np = parsin(2);
images.CSI.sr = images.CSI.sw / images.CSI.np; %spectral resolution, in Hz
images.CSI.srp = images.CSI.sr / ksp.ddf.transmit_freq; %spectral resolution, in ppm
images.CSI.soffp = ksp.ddf.ppm_ref; %spectral center, in ppm
% if isfield(images,'spbin') %adjust number of spectral points due to deletion
%     if isfield(images.spbin,'delete')
%         images.CSI.np = images.CSI.np - length(images.spbin.delete);
%     end
% end 
%Adjust some of the plotting parameters
images.plot.Cxbnds = [-1*images.CSI.fov(1),images.CSI.fov(1)] ...
    / 2 + images.CSI.off(1); %for 13C
images.plot.Cybnds = [-1*images.CSI.fov(2),images.CSI.fov(2)] ...
    / 2 + images.CSI.off(2);
images.plot.Czbnds = [-1*images.CSI.fov(3),images.CSI.fov(3)] ...
    / 2 + images.CSI.off(3);
end

% apodData: Adjusts apodization, apodizes 2D CSI data in spectral dimension
% and recalculates all images
%
function apodData(source,~)
apod = str2double(source.String);
recalcAll;
end

% gridShift: Shifts 2D CSI data to center nearest voxel on selected point
%
function gridShift(~,~)
set(ut,'String','Choose point to center nearest voxel on')
[xpos,ypos] = ginput(1); %select 1 point
% ID pH voxel on figure that contains selected point
xbnd = -.5*images.CSI.fov(1) + info.ppe + globshft(1); %original x- and  
    %y-edges of the FOV, prior to any shifting
ybnd = -.5*images.CSI.fov(2) + info.ppe2 + globshft(2);
for m = 1:images.CSI.npx(1)
    xl = xbnd + images.CSI.ps(1) * (m - 1);
    xu = xbnd + images.CSI.ps(1) * m;
    for n = 1:images.CSI.npx(2)
        yl = ybnd + images.CSI.ps(1) * (n - 1);
        yu = ybnd + images.CSI.ps(1) * n;
        if xpos > xl && xpos < xu && ypos > yl && ypos < yu
            dispvox = false(images.CSI.npx(1)); 
            dispvox(n,m) = true;
        end
    end
end
% Calculate distance between selected point and center of selected voxel
[yind,xind] = find(dispvox);
xposvox = xbnd + images.CSI.ps(1) * (xind - 0.5); %x-position 
    %of voxel center
yposvox = ybnd + images.CSI.ps(1) * (yind - 0.5); %y-position 
    %of voxel center
diff = [(xpos - xposvox),(ypos - yposvox)];
diffvox = diff ./ images.CSI.ps(1:2); %convert to # of 
    %voxels (voxel size after zerofilling)
%Reprocess CSI data with voxel shifts, plot
recalcAll;
end


%% INTERNAL FUNCTIONS - CALCULATION
%
% noiseDenoise: Adds noise to imaging data, performs tensor denoising,
% calculates maps, updates figure. Currently, performs a multi-denoise by
% apodizing + denoising steps, then undoing apodization and denoising again
%
function noiseDenoise(sizeLRin)
%Process CSI data with noise added (NO APODIZATION), perform 1st denoising
procCSI([images.CSI.npx(1) , images.CSI.np , 0 , diffvox(1) , ...
    diffvox(2)],'imgapodnoisy');
if nargin > 0
    sLR = sizeLRin;
else
    sLR = sizeLR;
end
sLRfirst = sLR;
sLRfirst(1) = 64/2;
[images.CSI.imgprogdenoised,nSVspec,nSVspat] = ...
    tensorDenoising(images.CSI.imgapodnoisy,sLRfirst);
% [images.CSI.imgprogdenoised,nSVspec,nSVspat] = ...
%     tensorDenoising(images.CSI.imgapodnoisy,sizeLR_firstpass);
dt = ksp.ddf.dwelltime/1000;
pad = (images.CSI.npx(1) - ksp.ddf.npix(1)) / 2;
t(1,1,:) = 0:dt:dt*(size(images.CSI.imgprogdenoised,3)-1); %time vector for 
    %apodization
% Apodize, then denoise data using tensor decomposition. Next, undo 
% apodization in order to return remaining noise to be approximately 
% Gaussian. Then, perform another denoising step with the same parameters.
temp = ifftn(images.CSI.imgprogdenoised); %to FID    
% apdmat = repmat(exp(t*apod/2*pi),ksp.ddf.npix(1), ...
%     ksp.ddf.npix(2),1); %exponential; an additional factor of 2*pi is 
%     %required to match SIVIC apodization results
sigma = 1 * sqrt(2 * log(2)) / pi / apod; %sigma parameter of 
    %Gaussian; an additional factor of 1/(2*pi) is required to match SIVIC 
    %apodization results
apdmat = repmat(exp(-1 * (t .^ 2 / 2 / (sigma ^ 2))),ksp.ddf.npix(1), ...
    ksp.ddf.npix(2),1); %gaussian
apdmat = padarray(apdmat,[1 1 0]*pad); %zerofill apdmat spatially 
temp = temp .* apdmat;
temp = fftn(temp); %back to spectra
[images.CSI.imgprogdenoised,nSVspec,nSVspat] = ...
    tensorDenoising(temp,sLR);
temp = ifftn(images.CSI.imgprogdenoised); %to FID
deapdmat = repmat(exp(1 * (t .^ 2 / 2 / (sigma ^ 2))),ksp.ddf.npix(1), ...
    ksp.ddf.npix(2),1); %gaussian     
deapdmat = padarray(deapdmat,[1 1 0]*pad); %zerofill deapdmat spatially    
temp = temp .* deapdmat;
temp = fftn(temp); %back to spectra
% images.CSI.imgdenoised1x = flip(images.CSI.imgdenoised1x,1);
[temp,nSVspec,nSVspat] = tensorDenoising(temp,sLR);
% Re-apodize data to match apodization of original dataset
temp = ifftn(temp); %to FID
% reapdmat = repmat(exp(-t*origapod*pi),ksp.ddf.npix(1), ...
%     ksp.ddf.npix(2),1); %exponential
sigma = 1 * sqrt(2 * log(2)) / pi / origapod; %sigma parameter of 
    %Gaussian
reapdmat = repmat(exp(-1 * (t .^ 2 / 2 / (sigma ^ 2))),ksp.ddf.npix(1), ...
    ksp.ddf.npix(2),1); %gaussian 
reapdmat = padarray(reapdmat,[1 1 0]*pad);
temp = temp .* reapdmat;
images.CSI.imgdenoisedFINAL = fftn(temp); %back to spectra
end

% calcMaps: Calculate integral, SNR, and linewidth maps. Input integer st 
% determines which maps are calculated (1 = all; 2 = all but original 
% image; 3 = all but original and apodized images; 4 = only noisy 
% original images and denoised image; 5 = denoised image only. If 
% deleteflg = true, spectral bins with deletions are used. If shiftcalcflg 
% = true, update shiftmap used for aligning spectra
% SPECTRAL DELETIONS AND PEAK ALIGNMENT NOT YET WORKING
%
function calcMaps(st)%,deleteflg,shiftcalcflg)
if nargin < 1
    st = 1;
end
% Calculate average noise standard deviation in all images. To correct for 
% any linear baseline slope, split noise region in half, reverse latter
% half, and add to former half. Then take std, multiply by sqrt(2) to
% account for noise reduction by summation
for ii = st:length(imgnames)
    nname = imgnames{ii};
    if ~contains(nname,'denoised') %don't procCSI for a denoised dataset
        if strcmp(nname,'img') || strcmp(nname,'imgnoisy') %use origapod
            procCSI([images.CSI.npx(1) , images.CSI.np , origapod , diffvox(1) , ...
                diffvox(2)],nname); 
        else %use apod
            procCSI([images.CSI.npx(1) , images.CSI.np , apod , diffvox(1) , ...
                diffvox(2)],nname);            
        end
    end
    halfpt = floor(length(spbin.noise)/2);
    sumnoise = images.CSI.(nname)(:,:,spbin.noise(1:halfpt)) + ...
        flip(images.CSI.(nname)(:,:,spbin.noise(halfpt+1:2*halfpt))...
        ,3);
    images.CSI.noisemap.(nname) = abs(std(sumnoise/2,0,3)) * sqrt(2);     
    % images.CSI.noisemap.(nname) = abs(std(...
    %     images.CSI.(nname)(:,:,spbin.noise),0,3));   
    % noise = reshape(images.CSI.(nname)(:,:,spbin.noise),1,[]);
    % images.CSI.noise.(nname) = abs(std(noise)); %takes from all voxels
    % Calculate maps of integral, SNR, linewidth over all images
    for jj = 1:length(spnames)
        spn = spnames{jj};
        maps.(nname).(spn).img = abs(...
            images.CSI.(nname)(:,:,spbin.(spn))); 
            %select spectral range to process, take abs val 
        maps.(nname).(spn).sum = sum(maps.(nname).(spn).img,3);
        maps.(nname).(spn).max = max(maps.(nname).(spn).img,[],3);
        maps.(nname).(spn).snr = maps.(nname).(spn).max ./ ...
            images.CSI.noisemap.(nname);
        %Calculate linewidths and midpoints of peaks (only for img and/or
        %imgapod)
        if strcmp(nname,'img') || strcmp(nname,'imgapod') 
            hm = maps.(nname).(spn).max / 2; %half-max value for each voxel
            for kk = 1:images.CSI.npx(1)
                for ll = 1:images.CSI.npx(2)
                    if maps.(nname).(spn).snr(ll,kk) > cutoff/2 %only calculate
                        %LW if SNR is high enough; use 1/2 the threshold for shift
                        %analysis
                        maxind = find(abs(images.CSI.(nname)(ll,kk,:)) == ...
                            maps.(nname).(spn).max(ll,kk)); %index of voxel max value
                        maxind = maxind(maxind >= spbin.(spn)(1) & ...
                            maxind <= spbin.(spn)(end)); %make sure 
                            %>1 value wasn't caught; make sure it's in spectral bin
                        hmind = zeros(2,1);
                        for mm = 1:round(3*length(spbin.(spn)))   
                            %move out from max index up and down, looking for half-max  
                            %values and corresponding indices
                            upind = maxind + mm;
                            dnind = maxind - mm;
                            if dnind < 1 %see if end was reached
                                if hmind(1) == 0
                                    hmind(1) = 1;
                                end
                            elseif abs(images.CSI.(nname)(ll,kk,dnind)) < hm(ll,kk)...
                                && hmind(1) == 0
                                %detect if peak has passed below half-max value for
                                %1st time
                                curdist = abs(abs(images.CSI.(nname)(ll,kk,dnind))...
                                    - hm(ll,kk)); %how close current point is to half-max
                                prevdist = abs(abs(images.CSI.(nname)(ll,kk,dnind+1))...
                                    - hm(ll,kk));%how close previous point is to half-max
                                hmind(1) = (curdist * (dnind + 1) + prevdist * dnind)...
                                    / (curdist + prevdist);
                                    %take mean value, weighting toward the index of
                                    %the point CLOSER to hm
                            end
                            if upind > size(images.CSI.(nname),3) %see if 
                                %end was reached
                                if hmind(2) == 0 
                                    hmind(2) = size(images.CSI.(nname),3);
                                end
                            elseif abs(images.CSI.(nname)(ll,kk,upind)) < hm(ll,kk)...
                                && hmind(2) == 0
                                %detect if peak has passed below half-max value for
                                %1st time
                                curdist = abs(abs(images.CSI.(nname)(ll,kk,upind))...
                                    - hm(ll,kk)); %how close current point is to half-max
                                prevdist = abs(abs(images.CSI.(nname)(ll,kk,upind-1))...
                                    - hm(ll,kk)); %how close previous point is to half-max
                                hmind(2) = (curdist * (upind - 1) + prevdist * upind)...
                                    / (curdist + prevdist);
                                    %take mean value, weighting toward the index of
                                    %the point CLOSER to hm
                            end                   
                            if sum(hmind > 0) == 2 %check if both half-max 
                                %indices have been found; break out of inner for() loop
                                break;
                            end
                        end  
                        maps.(nname).(spn).lw(ll,kk) = (hmind(2) - hmind(1)) ...
                            * images.CSI.sr; %peak linewidth, in Hz
                        maps.(nname).(spn).midpt(ll,kk) = mean(hmind); %midpoint 
                            %index inbetween half-max values (used for shift
                            %correction)
                    else
                        maps.(nname).(spn).lw(ll,kk) = 0;
                        maps.(nname).(spn).midpt(ll,kk) = 0;
                    end           
                end
            end
            if strcmp(nname,'img')
                mapname = ['Integral_' spn '_original'];
                maps.(mapname) = maps.(nname).(spn).sum; 
                lims.(mapname) = [min(reshape(maps.(mapname),1,[])) ...
                    max(reshape(maps.(mapname),1,[]))];
                mapname = ['LW_' spn '_original'];
                maps.(mapname) = maps.(nname).(spn).lw; 
                lims.(mapname) = [0 max(reshape(maps.(mapname) .* ...
                    (maps.(nname).(spn).snr > cutoff),1,[]))]; 
            elseif strcmp(nname,'imgapod')
                mapname = ['LW_' spn '_apodized'];
                maps.(mapname) = maps.(nname).(spn).lw; 
%                 lims.(mapname) = [0 max(reshape(maps.(mapname) .* ...
%                     (maps.(nname).(spn).snr > cutoff),1,[]))];
                lims.(mapname) = [0 median(reshape(...
                    maps.(mapname)(maps.(nname).(spn).snr > cutoff),1,[])) ...
                    + 2 * std(reshape(maps.(mapname)(maps.(nname).(spn).snr...
                    > cutoff),1,[]))];
            end
        end                
    end
    maps.(nname).pH = pKa + log10(maps.(nname).(spnames{1}).sum ./ ...
        maps.(nname).(spnames{2}).sum * fcorr); %pH map using magnitude integrals
%     maps.(nname).ratio = maps.(nname).(spnames{1}).max ./ ...
%         maps.(nname).(spnames{2}).max; %ratio map of peak heights  
end
maps.pH_original = maps.img.pH;
lims.pH_original = [6 8];
maps.pH_denoised = maps.imgdenoisedFINAL.pH;
lims.pH_denoised = lims.pH_original;
lims.(['LW_' spnames{1} '_original'])(2) = max(...
    [lims.(['LW_' spnames{1} '_original'])(2) lims.(['LW_' spnames{2} '_original'])(2)]);
lims.(['LW_' spnames{2} '_original']) = lims.(['LW_' spnames{1} '_original']);
lims.(['LW_' spnames{1} '_apodized'])(2) = max(...
    [lims.(['LW_' spnames{1} '_apodized'])(2) lims.(['LW_' spnames{2} '_apodized'])(2)]);
lims.(['LW_' spnames{2} '_apodized']) = lims.(['LW_' spnames{1} '_apodized']);
    %so both maps have same scaling
% Calculate change maps between original, noisy, and denoised
for ii = 1:length(spnames)
    spname = spnames{ii};
    mapname = ['SNRloss_' spname];
    maps.(mapname) = mean(reshape(maps.img.(spname).snr ./ ...
        maps.imgnoisy.(spname).snr,1,[]));
        %factor by which SNR decreased w/ noise addition
%     if ~nopeflg
%         disp([spname ': SNR decreased by factor of ' num2str(maps.(mapname),'%2.2f') ...
%             ' with noise addition'])
%     end
    mapname = ['UpSNR_' spname];
    maps.(mapname) = maps.imgdenoisedFINAL.(spname).snr ./ maps.imgnoisy.(spname).snr;% ...
%         .* (maps.imgdenoisedFINAL.(spname).snr >= cutoff);
        %factor by which SNR increased w/ denoising; as of 3/26/20 no longer
        %masking by SNR > cutoff threshold
    lims.(mapname) = [min(reshape(maps.(mapname),1,[])) max(reshape(maps.(mapname),1,[]))];
end
% maps.dRatio = maps.imgdenoisedFINAL.ratio ./ maps.img.ratio;% .* ...
% %     (maps.imgdenoisedFINAL.(spnames{1}).snr >= cutoff & maps.imgdenoisedFINAL.(spnames{2}).snr >= cutoff);  
%     % factor by which metabolite ratio changed; as of 4/17/20 no longer 
%     % masking by SNR > cutoff threshold in denoised dataset
% % lims.dRatio = [0 max(reshape(maps.dRatio,1,[]))];    
maps.dpH = maps.imgdenoisedFINAL.pH - maps.img.pH; %change in pH as result of  
    %ratio change; as of 4/17/20 no more SNR thresholding
% maps.dpH(isinf(maps.dpH)) = 0; %catch log10(0) voxels (due to SNR threshold)
% k = max(abs([min(reshape(maps.dpH,1,[])) max(reshape(maps.dpH,1,[]))]));
k = 0.3;
lims.dpH = [-k k]; %keep dpH lims centered on zero
maps.goodSNR = maps.imgdenoisedFINAL.(spnames{1}).snr >= cutoff & ...
    maps.imgdenoisedFINAL.(spnames{2}).snr >= cutoff; %voxels with good SNR, both peaks
lims.goodSNR = [0 1];
maps.goodpH = abs(maps.dpH) < tol; %voxels with |pH error| < tol in 
    %denoised compared with original data
lims.goodpH = [0 1];
maps.rescued_voxels = maps.goodSNR & maps.goodpH & ...
    (maps.imgnoisy.(spnames{1}).snr < cutoff | ...
    maps.imgnoisy.(spnames{2}).snr < cutoff );
    % indicates voxels in which SNR of one/both peaks increased above 
    % threshold as a result of apodizing + denoising, AND pH error < 
    % tolerance
lims.rescued_voxels = [0 1];  
maps.lost_voxels = ~maps.goodpH & (maps.imgnoisy.(spnames{1}).snr >= cutoff)...
    & (maps.imgnoisy.(spnames{2}).snr >= cutoff);
    % indicates voxels in which SNR of both peaks were good before and
    % after denoising, but pH is inaccurate after denoising
lims.lost_voxels = [0 1];
%New maps as of 3/27/20, to analyze pH and SNR in noisy data
maps.noisygoodSNR = maps.imgnoisy.(spnames{1}).snr >= cutoff & ...
    maps.imgnoisy.(spnames{2}).snr >= cutoff; %voxels with good SNR, both peaks,
    %in noisy data
% maps.noisydRatio = maps.imgnoisy.ratio ./ maps.imgapod.ratio;% .* ...
% %     maps.noisygoodSNR;  
%     % factor by which metabolite ratio changed; as of 4/17/20 no longer 
%     % masking based on SNR > cutoff in noisy dataset
maps.noisydpH = maps.imgnoisy.pH - maps.img.pH; %change in pH between 
    %noisy and original images; as of 4/17/20 no more SNR thresholding
maps.noisygoodpH = abs(maps.noisydpH) < tol; %voxels with |pH error| < tol in 
    %noisy compared with original data    
if ~nopeflg
    plotAxOvly;
end
% %Calculate map of spectral shifts required (in indices rather than ppm)
% %in order to align peaks, as follows. First, calculate midpoint 
% %inbetween peaks for each spectrum. Then take mean value and use to 
% %generate map of shifts for each spectra (to be used as input to
% %procCSI)
% if shiftcalcflg %calculate map for shifting spectra to align (used in procCSI)
%     midptmap = (maps.BiC.(nname).midpt + maps.CO2.(nname).midpt) ./ 2;
%     mpmean = mean(midptmap((maps.BiC.(nname).snr >= (cutoff/2)) & ...
%         (maps.CO2.(nname).snr >= (cutoff/2)))); %take mean of all midpoints 
%         %where both metabolites have SNR above 1/2 the cutoff
%     images.shiftmap = midptmap - mpmean; %map of shift values
%     images.shiftmap((maps.BiC.(nname).snr < (cutoff/2)) | ...
%         (maps.CO2.(nname).snr < (cutoff/2))) = 0; %don't shift spectra 
%         %with bad SNR, either metabolite
% end
end

% recalcAll: Using new inputs, noise and denoise, calculate all maps,
% redisplay spectra
%
function recalcAll(~,~)
    noiseDenoise; figure(100); calcMaps(1); plotSpecVox
end

% multCalc: Run noising+denoising with parameters multiple times, report
% statistics
%
function multCalc(sizeLRin,~)
% Initialize average-value maps
maps.avgdpH = zeros(images.CSI.npx(1));
maps.pctgoodSNR = zeros(images.CSI.npx(1)); %use to get % of runs of good   
    %SNR (BOTH peaks) per voxel
maps.pctgoodpHgoodSNR = zeros(images.CSI.npx(1)); %use to get % of runs of 
    %good pH per voxel GIVEN the voxel had sufficient SNR       
maps.pctrescue = zeros(images.CSI.npx(1)); %use to get % of runs per voxel
    %where denoising led to good pH + good SNR
maps.pctlost = zeros(images.CSI.npx(1)); %use to get % of runs per voxel
    %where denoising caused pH to be bad
% Initialize cumulative maps over all runs
for ii = 1:length(spnames)
    maps.(['avgUpSNR_' spnames{ii}]) = zeros(images.CSI.npx(1));
    maps.(['stdUpSNR_' spnames{ii}]) = zeros(images.CSI.npx(1));
    maps.cum.(['UpSNR_' spnames{ii}]) = zeros([images.CSI.npx(1) ...
        images.CSI.npx(1) maps.nRuns]);
    for jj = 1:length(imgnames)-2
        maps.(imgnames{jj+2}).(spnames{ii}).avgsnr = zeros(images.CSI.npx(1));
        maps.cum.(imgnames{jj+2}).(spnames{ii}).snr = zeros(...
            [images.CSI.npx(1) images.CSI.npx(1) maps.nRuns]);
    end
end
maps.cum.dpH = zeros([images.CSI.npx(1) images.CSI.npx(1) maps.nRuns]);
maps.cum.goodSNR = false([images.CSI.npx(1) images.CSI.npx(1) maps.nRuns]);
maps.cum.goodpHgoodSNR = false([images.CSI.npx(1) images.CSI.npx(1) ...
    maps.nRuns]);
%New maps as of 3/27/20, to analyze pH and SNR in noisy data
maps.cum.noisygoodSNR = false([images.CSI.npx(1) images.CSI.npx(1) maps.nRuns]);
maps.cum.noisygoodpHgoodSNR = false([images.CSI.npx(1) images.CSI.npx(1) ...
    maps.nRuns]);
nopeflg = true; %suppress text, figure plotting
% Conduct runs, store maps
for ii = 1:maps.nRuns
    if nargin > 0 && isnumeric(sizeLRin) %use specified value of sizeLR if included as input
        noiseDenoise(sizeLRin);
    else
        noiseDenoise; 
    end
    calcMaps(4); %maps for noisy + denoised images only
    for jj = 1:length(spnames)
        for kk = 1:length(imgnames)-3 %only imgnoisy and imgdenoisedFINAL
            maps.cum.(imgnames{kk+3}).(spnames{jj}).snr(:,:,ii) = ...
                maps.(imgnames{kk+3}).(spnames{jj}).snr;
        end
        maps.cum.(['UpSNR_' spnames{jj}])(:,:,ii) = maps.(['UpSNR_' spnames{jj}]);
    end
    maps.cum.dpH(:,:,ii) = maps.dpH;% .* maps.goodSNR; %as of 4/17/20 no 
        %longer thresholding by SNR
    maps.cum.goodSNR(:,:,ii) = maps.goodSNR;
    maps.cum.goodpHgoodSNR(:,:,ii) = maps.goodpH .* maps.goodSNR;
    maps.pctrescue = maps.pctrescue + maps.rescued_voxels;
    maps.pctlost = maps.pctlost + maps.lost_voxels;
    %New maps as of 3/27/20, to analyze pH and SNR in noisy data
    maps.cum.noisygoodSNR(:,:,ii) = maps.noisygoodSNR;
    maps.cum.noisygoodpHgoodSNR(:,:,ii) = maps.noisygoodpH .* maps.noisygoodSNR;
end
nopeflg = false;
for ii = 1:length(spnames)
    mapname = ['UpSNR_' spnames{ii}];
    mapname2 = ['avgUpSNR_' spnames{ii}];
    maps.(mapname2) = mean(maps.cum.(mapname),3);
%     lims.(mapname2) = [min(reshape(maps.(mapname2).*sneakmask,1,[])) ...
%         max(reshape(maps.(mapname2).*sneakmask,1,[]))];
    lims.(mapname2) = [min(reshape(maps.(mapname2),1,[])) ...
        max(reshape(maps.(mapname2),1,[]))];    
    if lims.(mapname2)(1) == lims.(mapname2)(2) %catch error: both equal
        lims.(mapname2)(2) = lims.(mapname2)(2) + 1;
    end
    for jj = 1:length(imgnames)-3 %only imgnoisy and imgdenoisedFINAL
        maps.(imgnames{jj+3}).(spnames{ii}).avgsnr = ...
            mean(maps.cum.(imgnames{jj+3}).(spnames{ii}).snr,3);
    end
    mapname3 = ['stdUpSNR_' spnames{ii}];
    maps.(mapname3) = std(maps.cum.(mapname),0,3);
end
maps.pctgoodSNR = mean(maps.cum.goodSNR,3) * 100;
maps.pctgoodSNRstd = std(maps.cum.goodSNR,0,3) * 100;
lims.pctgoodSNR = [0 100];
        
for ii = 1:images.CSI.npx(1) %for each voxel, we need to pull out ONLY the 
    %dpH values for which there is good SNR in order to calculate mean and 
    %stdev
    for jj = 1:images.CSI.npx(2)
        eval = maps.cum.dpH(jj,ii,:);
        pick = maps.cum.goodSNR(jj,ii,:);
        maps.avgdpH(jj,ii) = mean(eval(pick));
        maps.stddpH(jj,ii) = std(eval(pick));
    end
end
lims.avgdpH = [-0.2 0.2];
maps.pctgoodpHgoodSNR = mean(maps.cum.goodpHgoodSNR,3) * 100;% ...
%     ./ (maps.pctgoodSNR/100);
    %as of 4/17/20 includes ALL runs, not just ones when voxel had SNR > threshold
maps.pctgoodpHgoodSNR(isnan(maps.pctgoodpHgoodSNR)) = 0;
lims.pctgoodpHgoodSNR = [0 100];
maps.pctrescue = maps.pctrescue / maps.nRuns * 100;
lims.pctrescue = [0 100];
maps.pctlost = maps.pctlost / maps.nRuns * 100;
lims.pctlost = [0 100]; 
%New maps as of 3/27/20, to analyze pH and SNR in noisy data
maps.noisypctgoodSNR = mean(maps.cum.noisygoodSNR,3) * 100;
maps.noisypctgoodpHgoodSNR = mean(maps.cum.noisygoodpHgoodSNR,3) * 100;% ...
%     ./ (maps.noisypctgoodSNR/100);
    %as of 4/17/20 includes ALL runs, not just ones when voxel had SNR > threshold
maps.noisypctgoodpHgoodSNR(isnan(maps.noisypctgoodpHgoodSNR)) = 0;
plotAxOvly;
end

% svTest: Runs through different numbers of SVs kept in the denoising 
% (spectral, spatial) on multiply-selected voxels (NOT necessarily one 
% currently displayed!), plots max SNR increase (over all selected voxels)
% and % of good SNR + pH runs (over all selected voxels) as a function of 
% SVs. Runs through all combinations of SVs
%
function svTest(~,~)
brkflg = false; %so function doesn't break automatically, if ABORT button 
    %was clicked
set(ut,'String','SV test: Calculating...')
pause(.1)
[row,col] = find(selvox); %identify which voxel is currently selected
if isempty(row) || isempty(col) %ID if no voxel selected
    set(ut,'String','No voxel(s) chosen for test! Please select (a) voxel(s) and try again')
    brkflg = true;
else
    brkflg = false;    
end
% Iterate all combinations of spectral and spatial SVs, up until the specified
% nSVspecmax and nSVspatmax values
if ~brkflg
    %Disable GUI inputs that might mess up test
    set(snrs,'Enable','off');
    set(na,'Enable','off');
    set(snrs,'Enable','off');
    set(scsv,'Enable','off')
    set(stsv,'Enable','off')
    set(snrc,'Enable','off')
    set(nr,'Enable','off')
    set(ssvsm,'Enable','off')
    maxupSNR.avg = zeros(nSVspecmax,nSVspatmax);
%     maxupSNR.std = maxupSNR.avg; 
    pctgoodpHgoodSNR = maxupSNR.avg;
    for ii = 1:nSVspecmax
        for jj = 1:nSVspatmax
            drawnow; %update to see if user clicked ABORT button
            if brkflg
                brkflg = false; %reset for next time
                %Re-eable GUI inputs
                set(snrs,'Enable','on');
                set(na,'Enable','on');
                set(snrs,'Enable','on');
                set(scsv,'Enable','on')
                set(stsv,'Enable','on')
                set(snrc,'Enable','on')
                set(nr,'Enable','on')
                set(ssvsm,'Enable','on')
                set(ut,'String','SV test aborted.')
                return %exit function if ABORT button clicked
            end
            set(ut,'String',['SV test: Using ' num2str(ii) ' spectral SVs, '...
                num2str(jj) ' spatial SVs'])
            pause(.1)
            sLRin = [ii jj jj];
            multCalc(sLRin);
            % Find max SNR boost out of all voxels, spectral components
            pctvec = zeros(length(row),1);
            minSNRorig = pctvec;
            for kk = 1:length(row)
                SNRvec.(spnames{1})(kk) = maps.(['avgUpSNR_' spnames{1}])(row(kk),col(kk));
                SNRvec.(spnames{2})(kk) = maps.(['avgUpSNR_' spnames{2}])(row(kk),col(kk));
                SNRvec.max(kk) = max([SNRvec.(spnames{1})(kk),SNRvec.(spnames{2})(kk)]);
                    %SNR enhancement, larger of the two spectral components
                SNRvec.spmaxind(kk) = (SNRvec.(spnames{1})(kk) < ...
                    SNRvec.(spnames{2})(kk)) + 1; %index of spectral component 
                    %with a higher SNR enhancement
                pctvec(kk) = maps.pctgoodpHgoodSNR(row(kk),col(kk));
                minSNRorig(kk) = min([maps.imgnoisy.(spnames{1}).avgsnr(row(kk),col(kk)); ...
                    maps.imgnoisy.(spnames{2}).avgsnr(row(kk),col(kk))]);
            end
            maxind = find(SNRvec.max == max(SNRvec.max)); %index of voxel w/
                %max SNR boost
            maxpeakind = SNRvec.spmaxind(maxind); %index of spectral 
                %component with higher SNR boost
            maxupSNR.avg(ii,jj) = maps.(['avgUpSNR_' spnames{maxpeakind}])(row(maxind(1)),col(maxind(1)));
                %if a tie, just take the 1st one
%             maxupSNR.std(ii,jj) = maps.(['stdUpSNR_' spnames{maxpeakind}])(row(maxind(1)),col(maxind(1)));
            % Find minimum % of good pH + good SNR runs for all voxels
            minind = find(pctvec == min(pctvec)); %index of voxel w/
                %lowest % good pH + good SNR runs
            pctgoodpHgoodSNR(ii,jj) = maps.pctgoodpHgoodSNR(row(minind(1)),col(minind(1)));
                %if a tie, just take the 1st one
        end
    end
    % Find SV combinations that maximize both SNR boost and % of good pH
    % voxels
    [SNRrow,SNRcol] = find(maxupSNR.avg == max(reshape(maxupSNR.avg,1,[])));
    [pctrow,pctcol] = find(pctgoodpHgoodSNR == max(reshape(pctgoodpHgoodSNR,1,[])));
    if length(pctrow) > 1 %ID if multiple maxima were found
        multmaxflg = true;
        fltr = pctrow + pctcol;
        pctrow = pctrow(fltr == min(fltr)); %ID combination(s) that have lowest
            %total # of SVs kept
        pctcol = pctcol(fltr == min(fltr)); %ID combination(s) that have lowest
            %total # of SVs kept    
    else
        multmaxflg = false;
    end
    minSNRorig = min(minSNRorig);    
    % Plot figures in new window, display SV values that increase SNR boost
    % and % of good pH voxels
    figure; subplot(8,2,1:2); axis off
    title(['SV Test, noise amplitude = ' num2str(noiseamp,'%4.0f') ...
        ', min starting SNR = ' num2str(minSNRorig,'%2.1f') ', SNR cutoff = ' ...
        num2str(cutoff) ', ' num2str(maps.nRuns) ' runs'])
    subplot(8,2,[3 5 7 9 11]); surf(1:nSVspatmax,1:nSVspecmax,...
        maxupSNR.avg);
    xlabel('Spatial SVs kept'); ylabel('Spectral SVs kept');
    zlabel('Max SNR boost factor')    
    subplot(8,2,[4 6 8 10 12]); surf(1:nSVspatmax,1:nSVspecmax,...
        pctgoodpHgoodSNR);
    xlabel('Spatial SVs kept'); ylabel('Spectral SVs kept');
    zlabel('% of runs with good SNR + pH')
    subplot(8,2,[13 15]); plotAxOvly(true); axis square    
    subplot(8,2,16); axis off
    if multmaxflg %display only 1st discovered maximum, after filtering for 
        %lowest cumulative # of SVs
        title({['Max SNR boost with ' num2str(SNRrow) ' spectral SVs and ' ...
            num2str(SNRcol) ' spatial SVs'];['Max % of good pH with ' ...
            num2str(pctrow(1)) ' spectral SVs and ' num2str(pctcol(1)) ...
            ' spatial SVs']; '(additional maxima found)'});          
    else
        title({['Max SNR boost with ' num2str(SNRrow) ' spectral SVs and ' ...
            num2str(SNRcol) ' spatial SVs'];['Max % of good pH with ' ...
            num2str(pctrow) ' spectral SVs and ' num2str(pctcol) ' spatial SVs']});    
    end
    %Save data if specified
    if saveflg
        cd(save_dir)
        if isfolder(savename) == 0 %make directory if not already there in 
            %save path
            mkdir(savename)
        end
        cd(savename)
        basefilename = ['SVtest_' num2str(length(row)) 'voxels_zfill' ...
            num2str(spatzfillfactor) 'x_apod' num2str(apod,'%i') 'Hz_SNR' ...
            num2str(minSNRorig,'%2.1f') '_pHcutoff' num2str(tol,'%1.2f') ...
            '_SNRcutoff' num2str(cutoff) '_runs' num2str(maps.nRuns) ...
            '_upto' num2str(nSVspecmax) 'spec' num2str(nSVspatmax) 'spat'];
        ctr = 1;
        sname = basefilename;
        while 1 %this loop prevents saving over previous data
            if exist([sname '.mat'],'file')
                ctr = ctr + 1;
                sname = [basefilename '_' num2str(ctr,'%i')];
            else
                break;
            end
        end
        save([sname '.mat'],'maxupSNR','pctgoodpHgoodSNR','selvox')
        cd(home)
        set(ut,'String',['SV test: Results saved in ' save_dir '/' ...
            savename])
    else
        set(ut,'String','SV test complete! Figure may be behind this one')
    end
    %Re-eable GUI inputs
    set(snrs,'Enable','on');
    set(na,'Enable','on');
    set(snrs,'Enable','on');
    set(scsv,'Enable','on')
    set(stsv,'Enable','on')
    set(snrc,'Enable','on')
    set(nr,'Enable','on')
    set(ssvsm,'Enable','on')
    figure(100) %return to interactive figure for voxel selection
end
end

% snrTest: Runs through different noise values prior to denoising for 
% multiply-selected voxels, denoises using SV settings in GUI, plots max 
% SNR enhancement and min % of good SNR + pH runs (over all voxels and both
% spectral components) as a function of starting min SNR (over all voxels
% and both spectral components)
%
function snrTest(~,~)
set(ut,'String','SNR test: Calculating...')
pause(.1)
[row,col] = find(selvox); %identify which voxel(s) is (are) currently selected
if isempty(row) || isempty(col) %ID if no voxel selected
    set(ut,'String','No voxel chosen for test! Please select a voxel and try again')
    brkflg = true;
else
    brkflg = false;
end
% First determine amount of noise required to make lowest-SNR peak in all 
% voxels reach threshold, using apodized dataset. This will define the  
% middle value of noise values tested
if ~brkflg
    %Disable GUI inputs that might mess up test
    set(snrs,'Enable','off');
    set(na,'Enable','off');
    set(snrs,'Enable','off');
    set(scsv,'Enable','off')
    set(stsv,'Enable','off')
    set(snrc,'Enable','off')
    set(nr,'Enable','off')
    set(ssvsm,'Enable','off')    
    minSNRorig = zeros(length(row),1);
    pctvec = minSNRorig;
    pctnoisyvec = minSNRorig;
    s.String = num2str(cutoff);
    setSNR(s); %to get SNR to be measurable
    calcMaps(4); %to get image data using noiseamp for setting SNR
    for kk = 1:length(row)
%         minSNRorig(kk) = min([maps.img.(spnames{1}).snr(row(kk),col(kk)); ...
%             maps.img.(spnames{2}).snr(row(kk),col(kk))]);
        minSNRorig(kk) = min([maps.imgnoisy.(spnames{1}).snr(row(kk),col(kk)); ...
            maps.imgnoisy.(spnames{2}).snr(row(kk),col(kk))]);
    end
%     minind = find(minSNRorig == min(minSNRorig));
    minSNRorig = min(minSNRorig);
    noiseampSNR = noiseamp * minSNRorig / SNRset; %noise level required to 
        %set min SNR at desired SNR      
%     noiseampSNR = images.CSI.noisemap.img(row(minind),col(minind)) * ...
%         minSNRorig / cutoff; %noise level required to obtain min SNR at cutoff
    % Then generate vector of noise values
    nvec = 61;
    noisevec = logspace(log10(4*noiseampSNR),log10(0.2*noiseampSNR),nvec);
    % Run through different values of noise added to original image. Calculate
    % average min SNR over runs for each value of noise (to use with plotting
    % later)
    SNRvec = zeros(nvec,1);
    maxupSNR.avg = SNRvec;
    maxupSNR.std = SNRvec;
    pctgoodpHgoodSNR = SNRvec;
    noisypctgoodpHgoodSNR = SNRvec;
    nasave = noiseamp; %to save value so noiseamp gets returned to what it was 
        %at end
    for ii = 1:nvec
        drawnow; %update to see if user clicked ABORT button
        if brkflg
            brkflg = false; %reset for next time
            noiseamp = nasave; %return noiseamp to what it was
            %Re-eable GUI inputs
            set(snrs,'Enable','on');
            set(na,'Enable','on');
            set(snrs,'Enable','on');
            set(scsv,'Enable','on')
            set(stsv,'Enable','on')
            set(snrc,'Enable','on')
            set(nr,'Enable','on')
            set(ssvsm,'Enable','on')            
            set(ut,'String','SNR test aborted.')
            return %exit function if ABORT button clicked
        end
        set(ut,'String',['SNR test: Using noise value ' num2str(ii) ' of '...
            num2str(nvec)])
        pause(.1)    
        noiseamp = noisevec(ii);        
        multCalc;
        for kk = 1:length(row)
            minSNRvec.(spnames{1})(kk) = maps.imgnoisy.(spnames{1}).avgsnr(row(kk),col(kk));
            minSNRvec.(spnames{2})(kk) = maps.imgnoisy.(spnames{2}).avgsnr(row(kk),col(kk));
            minSNRvec.min(kk) = min([minSNRvec.(spnames{1})(kk),minSNRvec.(spnames{2})(kk)]);
                %starting SNR, smaller of the two spectral components
            minSNRvec.spminind(kk) = (minSNRvec.(spnames{1})(kk) > ...
                minSNRvec.(spnames{2})(kk)) + 1; %index of spectral component 
                %with a lower starting SNR   
            pctvec(kk) = maps.pctgoodpHgoodSNR(row(kk),col(kk));
            pctnoisyvec(kk) = maps.noisypctgoodpHgoodSNR(row(kk),col(kk));
        end   
        minind = find(minSNRvec.min == min(minSNRvec.min)); %index of voxel
            %with lowest starting SNR
        minpeakind = minSNRvec.spminind(minind); %index of spectral 
                %component with lowest starting SNR    
        SNRvec(ii) = maps.imgnoisy.(spnames{minpeakind}).avgsnr(row(minind),col(minind));       
        maxupSNR.avg(ii) = maps.(['avgUpSNR_' spnames{minpeakind}])(row(minind),col(minind));
        maxupSNR.std(ii) = maps.(['stdUpSNR_' spnames{minpeakind}])(row(minind),col(minind));
        pctgoodpHgoodSNR(ii) = min(pctvec);
        noisypctgoodpHgoodSNR(ii) = min(pctnoisyvec);
    end
    noiseamp = nasave; %return noiseamp to what it was
    % Plot figures in new window
%     figure; subplot(10,2,1:2); axis off
    figure; subplot(7,2,1:2); axis off
    title(['SNR Test, ' num2str(sizeLR(1)) ' spectral SVs, ' ...
        num2str(sizeLR(2)) ' spatial SVs, SNR cutoff = ' num2str(cutoff) ...
        ', ' num2str(maps.nRuns) ' runs'])
    stdplot = [(maxupSNR.avg - maxupSNR.std),...
        (2 * maxupSNR.std)];    
    subplot(7,2,[3 5 7]); 
    q = area(SNRvec,stdplot,'LineStyle',':','FaceAlpha',0.2);
    q(1).FaceAlpha = 0;
    q(1).EdgeAlpha = 0;
    q(2).FaceColor = 'g'; 
    hold on; plot(SNRvec,maxupSNR.avg); %plotting this second allows user 
        %to better ID average values of SNR increase from plot
    xlabel('Min SNR prior to denoising'); ylabel('Max SNR boost factor')
%     subplot(10,2,[10 12 14]); plot(SNRvec,noisypctgoodpHgoodSNR)
    subplot(7,2,[4 6 8]); plot(SNRvec,noisypctgoodpHgoodSNR)
    hold on; plot(SNRvec,pctgoodpHgoodSNR);
    xlabel('Min SNR prior to denoising'); ylabel('% of runs with good SNR + pH')
    legend('Before denoising','After denoising');  
%     subplot(10,2,[3 5 7]); plot(noisevec,maxupSNR.avg); hold on;
%     p = area(noisevec,stdplot,'LineStyle',':','FaceAlpha',0.2);
%     p(1).FaceAlpha = 0;
%     p(1).EdgeAlpha = 0;
%     p(2).FaceColor = 'g'; 
%     xlabel('Noise amplitude added'); set(gca,'xdir','reverse'); 
%     ylabel('Max SNR boost factor')
%     subplot(10,2,[4 6 8]); plot(noisevec,noisypctgoodpHgoodSNR)
%     hold on; plot(noisevec,pctgoodpHgoodSNR);
%     xlabel('Noise amplitude added'); set(gca,'xdir','reverse');
%     ylabel('% of runs with good SNR + pH'); 
%     legend('Before denoising','After denoising');
%     subplot(10,2,[9 11 13]); plot(SNRvec,maxupSNR.avg); hold on;    
%     subplot(10,2,[17 19]); plotAxOvly(true); axis square
    subplot(7,2,[11 13]); plotAxOvly(true); axis square
    %Save data if specified
    if saveflg
        cd(save_dir)
        if isfolder(savename) == 0 %make directory if not already there in 
            %save path
            mkdir(savename)
        end
        cd(savename)
        basefilename = ['SNRtest_' num2str(length(row)) 'voxels_zfill' ...
            num2str(spatzfillfactor) 'x_apod' num2str(apod,'%i') ...
            'Hz_pHcutoff' num2str(tol,'%1.2f') '_SNRcutoff' num2str(cutoff) ...
            '_runs' num2str(maps.nRuns) '_' num2str(sizeLR(1)) 'specSVs_' ...
            num2str(sizeLR(2)) 'spatSVs'];
        ctr = 1;
        sname = basefilename;
        while 1 %this loop prevents saving over previous data
            if exist([sname '.mat'],'file')
                ctr = ctr + 1;
                sname = [basefilename '_' num2str(ctr,'%i')];
            else
                break;
            end
        end
        save([sname '.mat'],'maxupSNR','SNRvec','noisypctgoodpHgoodSNR',...
            'pctgoodpHgoodSNR','selvox')
        cd(home)
        set(ut,'String',['SNR test: Results saved in ' save_dir '/' ...
            savename])
    else
        set(ut,'String','SNR test complete! Figure may be behind this one')
    end
    %Re-eable GUI inputs
    set(snrs,'Enable','on');
    set(na,'Enable','on');
    set(snrs,'Enable','on');
    set(scsv,'Enable','on')
    set(stsv,'Enable','on')
    set(snrc,'Enable','on')
    set(nr,'Enable','on')
    set(ssvsm,'Enable','on')
    figure(100) %return to interactive figure for voxel selection
end
end


%% INTERNAL FUNCTIONS - SET VALUES
%

% setNoise: Sets amplitude of noise standard deviation, adds to original
% CSI data
%
function setNoise(source,~)
    noiseamp = str2double(source.String);
    set(ut,'String','')
    recalcAll;
end

% setSNR: Calculates and sets amplitude of noise standard deviation in 
% order to obtain specified SNR for lowest-SNR metabolite in DISPLAYED 
% voxel (apodized spectrum)
%
function setSNR(source,~)
[row,col] = find(dispvox); %identify which voxel is currently selected for display
if isempty(row) || isempty(col) %ID if no voxel selected
    set(ut,'String','No voxel chosen for test! Please select a voxel and try again')
    brkflg = true;
else
    brkflg = false;    
end
if ~brkflg    
    SNRset = str2double(source.String);
    minSNRorig = min([maps.img.(spnames{1}).snr(row,col) ...
        maps.img.(spnames{2}).snr(row,col)]);
    noiseamp = images.CSI.noisemap.img(row,col) * minSNRorig / SNRset; %noise 
        %level required to obtain min SNR at desired SNR    
%     minSNRorig = min([maps.imgnoisy.(spnames{1}).snr(row,col) ...
%         maps.imgnoisy.(spnames{2}).snr(row,col)]);
%     noiseamp = noiseamp * minSNRorig / SNRset; %noise level required to 
%         %set min SNR at desired SNR    
    recalcAll;
    set(ut,'String',['Voxel [x = ' num2str(col) ', y = ' num2str(row) ...
        '] from top lefthand corner: noise set for desired SNR']) 
    set(na,'String',num2str(noiseamp))
end
end

% setSpecSV: Sets number of SPECTRAL singular values to keep during
% denoising
%
function setSpecSV(source,~)
    sizeLR(1) = str2double(source.String);
end

% setSpatSV: Sets number of SPATIAL singular values to keep during
% denoising (equal for both dimensions)
%
function setSpatSV(source,~)
    sizeLR(2) = str2double(source.String);
    sizeLR(3) = sizeLR(2);
end

% setCutoff: Sets SNR cutoff for thresholding maps
%
function setCutoff(source,~)
    cutoff = str2double(source.String);
end

% setRuns: Sets number of runs for multirun
%
function setRuns(source,~)
    maps.nRuns = str2double(source.String);
    set(mr,'String',['MULTIRUN (x' num2str(maps.nRuns,'%i') ')']); %update 
        %button label
end

% setSVspecmax: Sets number of SVs (spectral) to run for SV test
%
function setSVspecmax(source,~)
    nSVspecmax = str2double(source.String);
end

% abortTest: Break out of currently-running SV or SNR test
%
function abortTest(~,~)
brkflg = true;
end

% dispVoxOnOff: Toggle voxel selection on and off
%
function dispVoxOnOff(~,~)
    finishflg = ~finishflg;
    if finishflg
        set(tb,'String','Turn individual voxel display on')
    else
        set(tb,'String','Turn individual voxel display off')
    end
    dispVox;
end

% axSelect: Plot selected axial slice
%
function axSelect(source,~)
axsl = round(source.Value);
plotAxOvly;
end

% selOvly: Set active overlay for display purposes
%
function selOvly(source,~)
ovlyname = ovlynames{source.Value};
if contains(ovlyname,'LW_') %display linewidth results for selected voxels
    if contains(ovlyname,'_original')
        nname = 'img';
    elseif contains(ovlyname,'_apodized')
        nname = 'imgapod';
    end
    [row,col] = find(selvox); %identify which voxel(s) is (are) currently selected
    if ~isempty(row) && ~isempty(col) %verify that voxels were selected
        for ii = 1:length(spnames)
            spname = spnames{ii};
            for jj = 1:length(row)
                LWs.(spname).vals(jj) = maps.(nname).(spname).lw(row(jj),col(jj));
            end
            LWs.(spname).avg = mean(LWs.(spname).vals);
            LWs.(spname).diff = max(LWs.(spname).vals) - min(LWs.(spname).vals);
        end
        if LWs.(spnames{1}).diff > LWs.(spnames{2}).diff %ID which has the 
                %greater difference in linewidth
            maxpctLW = LWs.(spnames{1}).diff / LWs.(spnames{1}).avg * 100; 
                %normalize by mean linewidth across all voxels
        else
            maxpctLW = LWs.(spnames{2}).diff / LWs.(spnames{2}).avg * 100;
        end
        set(ut,'String',[nname ': max % difference in linewidth, relative' ...
            ' to mean across voxels = ' num2str(maxpctLW,'%2.0f') '%'])
    end
end
plotAxOvly;
end


% setOvlyTransparency: Adjust transparency of all 13C overlay
%
function setOvlyTransparency(source,~)
otp = source.Value;
plotAxOvly;
end

% setZoom: Zoom in on all images in final figure
%
function setZoom(source,~)
zoomfac = source.Value;
plotAxOvly;
end

% resetZoom: Reset zoom to original amount
%
function resetZoom(~,~)
zoomfac = 0;
plotAxOvly;
end

% setMultVoxNum: Allows user to set number of voxels to select
%
function setMultVoxNum(source,~)
nvox = str2double(source.String);
end

% dispMultVoxOnOff: Toggle outlines of multiple selected voxels on and off
%
function dispMultVoxOnOff(~,~)
    multdispflg = ~multdispflg;
    if multdispflg
        set(tbmv,'String','Turn voxel outlines off')
    else
        set(tbmv,'String','Turn voxel outlines on')
    end
    plotAxOvly;
end

% setSave: Toggle whether SV/SNR test results are saved to .mat files
%
function setSave(source,~)
if source.Value == 1
    saveflg = true;
else
    saveflg = false;
end
end

% finishFig: Finish all figure stuff, continue function
%
function finishFig(~,~)
finishflg = true;
close all
end


%% INTERNAL FUNCTIONS - IMAGE DISPLAY, INTERACTION
%
% plotAxOvly: plots axial 1H image + overlay voxels
%
function plotAxOvly(lessflg)
if nargin < 1 
    lessflg = false;
end
if gcf == voxfig || lessflg %only plot if voxfig is currently selected, or 
    %if displaying voxels on SV/SNR test figs
    % First, convert image to grayscale (necessary for overlays)
    plotimg = squeeze(images.ref.img(:,:,axsl));
    hlims = [min(min(plotimg)) max(max(plotimg))];
    plotimg = repmat(mat2gray(plotimg,hlims),[1 1 3]);
    if ~lessflg
        subplot(1,2,1); 
    end
    imagesc(images.plot.Hxbnds,images.plot.Hybnds,plotimg);
    hold on
    if ~lessflg
        % Plot overlay
        v = imagesc((images.plot.Cxbnds-images.CSI.off(1))*pixrat+images.CSI.off(1),...
            (images.plot.Cybnds-images.CSI.off(2))*pixrat+images.CSI.off(2),...
            maps.(ovlyname),lims.(ovlyname));
        colormap(jet); colorbar; axis('square');
        % set(v,'AlphaData',ones(images.CSI.npx(1))*otp.*sneakmask);
        set(v,'AlphaData',ones(images.CSI.npx(1))*otp);
    end
    % Outline any voxels selected for SV/SNR tests
    if multdispflg || lessflg
        [olvoxy,olvoxx] = find(selvox);
        for ii = 1:length(olvoxx) %plot lines around each selected voxel
            xl = images.plot.Cxbnds(1) + images.CSI.ps(1) * (olvoxx(ii)-1);
            xu = images.plot.Cxbnds(1) + images.CSI.ps(1) * olvoxx(ii);
            yl = images.plot.Cybnds(1) + images.CSI.ps(1) * (olvoxy(ii) - 1);
            yu = images.plot.Cybnds(1) + images.CSI.ps(1) * olvoxy(ii);
            if lessflg
                plot([xl xu xu xl xl],[yl yl yu yu yl],'g-','LineWidth',1);
            else
                plot([xl xu xu xl xl],[yl yl yu yu yl],'g-','LineWidth',3);
            end
        end
    end
    % Plot and label scale bar
    plotaxes = [images.plot.Hxbnds,images.plot.Hybnds] .* ones(1,4) ...
        * (1-zoomfac);
    if ~lessflg %don't plot scale bar if miniplot
        barx = [plotaxes(2)*0.9-4 plotaxes(2)*0.9];
        % bary = [plotaxes(4) plotaxes(4)] * 0.9;
        bary = [14 14];
        plot(barx,bary,'w-','LineWidth',4)
        text(mean(barx),bary(1)*.8/.9,'4 mm','FontSize',16,'Color','w',...
             'HorizontalAlignment','center');
    end
    axis(plotaxes); %adjust zoom
    hold off
end
end

% plotSpecVox: Plots CSI spectra from selected voxel
%
function plotSpecVox
delete(spax); %clear currently plotted spectra
clear spax
if ~isempty(dispvox(dispvox)) %check to see if any voxels selected for plotting
    [sprow,spcol] = find(dispvox); %get coordinates for spectra to plot
    xl = .25;
    yl = xl;
    xp = .5 + [0 xl/1.7 xl*2/1.7];
    yp = .9 - [yl yl*2.5];
    %plot original spectrum first
    spax(1) = axes('Position',[xp(1) yp(1) xl yl]);
    plot(abs(squeeze(images.CSI.img(sprow,spcol,:))));
    title('img'); axis('square'); axis([1 images.CSI.np 0 max(abs(...
        reshape(images.CSI.img,1,[]))) * 1.1]);
    set(gca,'Xtick',[]); set(gca,'Ytick',[]);     
    spax(2) = axes('Position',[xp(2) yp(1) xl yl]);
    plot(abs(squeeze(images.CSI.imgnoisy(sprow,spcol,:))));
    title('imgnoisy'); axis('square'); axis([1 images.CSI.np 0 max(abs(...
        reshape(images.CSI.imgnoisy,1,[]))) * 1.1]);
    set(gca,'Xtick',[]); set(gca,'Ytick',[]);    
    spax(3) = axes('Position',[xp(3) yp(1) xl yl]);
    plot(abs(squeeze(images.CSI.imgapod(sprow,spcol,:))));
    title('imgapod'); axis('square'); axis([1 images.CSI.np 0 max(abs(...
        reshape(images.CSI.imgapod,1,[]))) * 1.1]);
    set(gca,'Xtick',[]); set(gca,'Ytick',[]);       
    spax(4) = axes('Position',[xp(1) yp(2) xl yl]);
    plot(abs(squeeze(images.CSI.imgapodnoisy(sprow,spcol,:)))); 
    title('imgapodnoisy'); axis('square'); axis([1 images.CSI.np 0 max(abs(...
        reshape(images.CSI.imgapodnoisy,1,[]))) * 1.1]);
    set(gca,'Xtick',[]); set(gca,'Ytick',[]);
    spax(5) = axes('Position',[xp(2) yp(2) xl yl]);
    plot(abs(squeeze(images.CSI.imgprogdenoised(sprow,spcol,:)))); 
    title('imgprogdenoised'); axis('square'); axis([1 images.CSI.np 0 max(abs(...
        reshape(images.CSI.imgprogdenoised,1,[]))) * 1.1]);
    set(gca,'Xtick',[]); set(gca,'Ytick',[]);    
    spax(6) = axes('Position',[xp(3) yp(2) xl yl]);
    plot(abs(squeeze(images.CSI.imgdenoisedFINAL(sprow,spcol,:)))); 
    title('imgdenoisedFINAL'); axis('square'); axis([1 images.CSI.np 0 max(abs(...
        reshape(images.CSI.imgdenoisedFINAL,1,[]))) * 1.1]);
    set(gca,'Xtick',[]); set(gca,'Ytick',[]);    
%         for m = 1:2
%             text(images.CSI.np*1/4,lims.plot(2)*(1.1-.065*(2+m)),[spnames{m} ...
%                 ': Avg SNR = ' num2str(maps.(imgnames{l+1}).(spnames{m}).avgsnr(sprow,spcol),'%3.3f')]);
%         end
else
    spax = [];
end
end

% dispVox: Allows user to select voxel for seeing spectra
%
function dispVox(~,~)
while 1
    drawnow; %update variables
    if ~finishflg || isvalid(voxfig) %this loop breaks when finishflg is 
        %true or if GUI figure has been deleted    
        try
            [xpos,ypos] = ginput(1); %select 1 voxel
        catch
        end
    end
    drawnow; %update variables
    if finishflg || ~isvalid(voxfig) %this loop breaks when finishflg is 
        %true or if GUI figure has been deleted
        return;
    else
        % ID selected pH voxels on figures
        for m = 1:images.CSI.npx(1)
            xl = images.plot.Cxbnds(1) + images.CSI.ps(1) * (m - 1);
            xu = images.plot.Cxbnds(1) + images.CSI.ps(1) * m;
            for n = 1:images.CSI.npx(2)
                yl = images.plot.Cybnds(1) + images.CSI.ps(1) * (n - 1);
                yu = images.plot.Cybnds(1) + images.CSI.ps(1) * n;
                if xpos > xl && xpos < xu && ypos > yl && ypos < yu
                    dispvox = false(images.CSI.npx(1)); 
                    dispvox(n,m) = true;
                end
            end
        end
    % Plot selected CSI voxel on the other half of the figure
    plotSpecVox;
    end
end
end

% selMultVox: Allows user to select multiple voxel for SV and SNR tests
%
function selMultVox(~,~)
    selvox = false(images.CSI.npx(1)); %reset selvox
    multdispflg = true; %so selected voxels can be seen during selection
    set(tbmv,'String','Turn voxel outlines off')
    plotAxOvly; %remove voxel outlines from last time
    % ID selected pH voxels on figures
    for p = 1:nvox
        set(ut,'String',['Please select ' num2str(nvox-p+1,'%i') ' voxel(s)']);
        [xpos,ypos] = ginput(1); %select 1 voxel, nvox times
        for m = 1:images.CSI.npx(1)
            xl = images.plot.Cxbnds(1) + images.CSI.ps(1) * (m - 1);
            xu = images.plot.Cxbnds(1) + images.CSI.ps(1) * m;
            for n = 1:images.CSI.npx(2)
                yl = images.plot.Cybnds(1) + images.CSI.ps(1) * (n - 1);
                yu = images.plot.Cybnds(1) + images.CSI.ps(1) * n;
                if xpos > xl && xpos < xu && ypos > yl && ypos < yu
                    selvox(n,m) = true;
                end
            end
        end
        plotAxOvly; %update figure to show selected voxel
    end
    set(ut,'String','');
end


end


%% EXTERNAL FUNCTIONS
%
% refimgorder: Reorders 1H FIDs from raw data
function ksp = refimgorder(re,im,np,etl)
ksp = re + 1i * im;
if nargin < 4 %fsems needs echo train length as input; otherwise, etl = 1
    etl = 1;
end
pelist = zeros(etl*2,1);
for ii = 1:etl
    pelist(ii) = 2 * (etl - ii) + 1;
    pelist(etl+ii) = 2 * ii;
end
ksp = reshape(ksp,[np(2) etl np(3) np(1)/etl]);
ksp = permute(ksp,[1 4 2 3]);
ksp = reshape(ksp,[np(2) np(1)/etl/2 etl*2 np(3)]);
ksp = ksp(:,:,pelist,:);
ksp = reshape(ksp,[np(2) np(1) np(3)]);
ksp = flip(ksp,1); %need to flip images along x and y
ksp = flip(ksp,2);
end

% tensorDenoising - Written by D. Korenchan, based upon MATLAB code by 
% C. Baligand, C. Najac, and J. Brender/Galen Reed
% --    Utilizes tensor value decomposition methods, as described in: 
%       Brender et al, BioRxiv 2018
%
% Reads in data file (using SIVIC systems) and performs tensor
% denoising. Imaging data are saved in DDF format
function [imgout,nSVspec,nSVspat] = tensorDenoising(image,sizeLR)
% Perform tensor value decomposition upon data
%
% dataSize = size(image);
% options for phase correction search (if used)
% bruteForce = 0;
% simplexZeroOrder = 1;
% simplexZeroAndFirstOrder = 2;

% params.phaseSensitive = false;
% params.phaseCorrectTimePoints = false; 
% params.nonNegativePenalty = true;
% params.searchMethod = simplexZeroOrder;
% params.baselineCorrect = false;
% params.specPoints = dataSize(1);
% params.timePoints = dataSize(2);
% params.channels = dataSize(3);

spectra = permute(image,[3 1 2]); %need spectral dimension to be 1st
% spectra = fftshift(fftn(fftshift(images.img)));
% phaseCorrectedSpectra = phaseAndBaselineCorrect(spectra, params);


[U,S,sv] = mlsvd(spectra);

% Plot singular values obtained from tensor decomposition
%
nSVspec = length(sv{1});
nSVspat = length(sv{2});
% figure();
for ii = 1:length(sv)
%     y = sv{ii};
%     subplot(1,length(sv),ii);
%     semilogy(y);
%     title(['Singular Values, Dimension ' num2str(ii)])
% Truncate matrices, reconstruct data
%
%     disp(['Number of singular values kept, dimension ' num2str(ii) ': ' num2str(sizeLR(ii))]);
    Utrunc{ii} = U{ii}(:, 1:sizeLR(ii));
end
if length(sv) == 2
    Strunc = S(1:sizeLR(1), 1:sizeLR(2));
elseif length(sv) == 3
    Strunc = S(1:sizeLR(1), 1:sizeLR(2), 1:sizeLR(3));
elseif length(sv) == 4
    Strunc = S(1:sizeLR(1), 1:sizeLR(2), 1:sizeLR(3), 1:sizeLR(4));
elseif length(sv) == 5
    Strunc = S(1:sizeLR(1), 1:sizeLR(2), 1:sizeLR(3), 1:sizeLR(4), 1:sizeLR(5));
end
specsLR = lmlragen(Utrunc, Strunc);
imgout = permute(specsLR,[2 3 1]); %get back to y by x by f
% imgout = ifftshift(ifftn(ifftshift(specsLR))); %return to k-space data
end