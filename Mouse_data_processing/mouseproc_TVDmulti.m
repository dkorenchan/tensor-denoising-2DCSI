function mouseproc_TVDmulti(imageprev,roiprev)
% mouseproc.m - Written by D. Korenchan, based upon ASL3.m by C. Baligand
%
% Processes and displays 1H axial + coronal images and HP 13C bicarbonate  
% 2D CSI data, along with associated overlays, using a denoising algorithm 
% based upon the Tucker decomposition. The user is asked to specify files 
% to load as well as processing parameters, including zerofilling, 
% spectral apodization, grid shifting, integral bounds for peak area 
% calculation and SNR thresholding (interactive display is used to identify
% best spectral apodization settings). The user selects voxels on the 1H 
% axial and coronal images for any number of ROIs, and the script will 
% calculate and return the sum, mean, standard deviation, maximum, and 
% minimum from each ROI. Images and ROI data may be saved at the end of the 
% script as MATLAB variables and/or text file output.
%
%
% INPUTS (BOTH OPTIONAL)
% prev_images:  Structure containing all 1H and 13C images and
%               processing/acquisition parameters from running 
%               mouseproc_TVD.m before
% prev_roi:     Structure containing ROI geometric parameters and
%               statistics, previously processed
%
% OUTPUTS
% images:  Structure containing 1H axial + coronal and 13C bicarbonate 
%          images, along with processing/acquisition parameters         
% rois:    Structure containing ROI geometric parameters and
%          statistics      
%
% USAGE - DATA PROCESSING (NO INPUTS):
% 1.    When prompted, specify .fid of the 2D CSI data (script may also ask
%       for specification of axial + coronal 1H image files)
% 2.    When prompted, specify zerofilling, apodization, integral bounds, 
%       and grid shifting parameters
% 3.    When prompted, specify noise region of spectra to obtain noise
%       standard deviation
% 4.    When interactive GUI displays, use to specify ROI names and select
%       voxels
% 5.    When prompted, specify folder and .mat filename for saving
%       variables, or click Cancel to skip saving
%
% CHILD SCRIPTS:    load_fids.m, load_procpar.m,
%                   scripts from SIVIC software
% UPDATES:
%   4/29/20:    Ported over from mouseproc.m
%   5/4/20:     Got function to go all the way to end. Haven't added in GUI
%               feature at end to compare spectra. 
%   5/5/20:     Added ability to specify different integral bounds for
%               BiC/CO2 between original and denoised datasets, as well as
%               the ability to visualize these bounds. 
%   5/6/20:     Fixed image loading (ROI loading already worked); also
%               improved display of integral bounds and completed spectral
%               voxel display in final figure
%   5/6/20:     Included additional processing option: user is able to cut
%               out any number of spectral regions that may contain 
%               nuisance peaks prior to apodization. Function cuts them
%               out, inverse FTs, apodizes, and FTs again to regenerate the
%               spectra.
%   5/8/20:     Completed debugging of spectral eliminations, FT/inverse
%               FT, etc. Verified that spectra are nearly identical to
%               SIVIC. 
%   5/8/20:     Fixed linewidth calculation bug: script was unable to find
%               half-max points outside of the defined spectral bin
%   5/10/20:    Finished fixes to linewidth calculation. Also fixed noise
%               bin problem: if deletions occurred in spectrum, noise bin
%               could shift over and ruin SNR calculations
%   5/11/20:    Added spectral alignment option to shift all spectral FIDs
%               so as to align the peaks
%   5/12/20:    Made spectral deletions defined on spectra WITHOUT
%               apodization
%   5/24/20:    BIG BUG FIX: Spatial grid shifting has not been being 
%               performed since code was changed to do k-space processing 
%               independent of SIVIC!! Added code to shift data in k-space 
%               based upon diffvox variable
%   5/26/20:    Added triple-Tucker decomposition, and moved CSI data
%               loading to be outside procCSI subfunction. 
%   5/27/20:    MAJOR BUG FIX: Spectral shifting factor tensor needs to be
%               FT'd spatially (both dimensions) to work properly! Also
%               found that there was some intrinsic pH bias due to
%               differences in processing b/w procCSI and operations in the
%               main function inbetween denoising steps - moved all
%               denoising to procCSI and verified that minimal pH bias 
%               (<0.001 pH unit) is introduced by processing operations
%   5/28/20:    Shifting must be done prior to denoising, which is done
%               prior to zerofilling, so peak shifting parameters are
%               calculated on native resolution data (which are also
%               displayed in first GUI). Also switched back to noise
%               calculation w/o addition, and made peak alignment based on
%               BiC only
%   6/3/20:     Troubleshooting for peak alignment + proper matching with
%               SIVIC: Discovered that SIVIC shifts spectra over ~1/2
%               voxel, both x and y! Corrected for this - now they match
%               pretty darn well!! Everything seems to check out now 
%               (though still not convinced the peak alignment is working
%               perfectly...maybe due to inaccuracies in ID'ing peak middle)
%   6/3/20:     Removed first denoise; found that even to 3/4 spectral SVs
%               that pH bias as high as 0.05 unit introduced! (MS494)
%   6/6/20:     Tried progressive-denoise Tucker in 25 Hz increments, since
%               double-/triple-Tucker weren't working well
%

%% VALUES TO ADJUST PRIOR TO RUNNING
%
BiCppm = [210.8,209]; %integral bounds for BiC, to use with non-denoised 
    %data
CO2ppm = [175.3,172.9]; %integral bounds for CO2, to use with non-denoised 
    %data
images.CSI.origapod = 10; %apodization to use with non-denoised data, in Hz
images.CSI.apodincr = 25; %max increment (in Hz) to apodize each time
sizeLR = [5, 8, 8]; %core tensor lengths to KEEP for each dimension 
    %(spectral, spatial, spatial) when denoising
load_path = '/Users/sf865719/Lab/Data/600_Data/pH_invivo/GLC/TRAMP'; %path to data
save_dir = '/Users/sf865719/Lab/MATLAB/Image_processing/Saved_data/mouseproc';
    %path for saving data
svk_path = '/Applications/SIVIC.app/Contents/sivic/local/bin/'; %path to 
    %SIVIC scripts
brukflg = false; %set to true if loading from Bruker ser file
filestrref = {'sems_mouse_axial','sems_mouse_cor'}; %subset of reference image filename
% filestrref = {'ref','cor'}; %subset of reference image filename
lastflg = true; %selects final file found when searching with filestrref;
    %otherwise, specified by user
globshft = [0,0]; %[mm right,mm down] shift of entire 13C data relative to 1H 
    %(axial plane)   
images.pHthr = 3; %SNR factor used for thresholding, 13C pH
images.fcorr = 9.276; %flip angle correction factor for pH imaging


%% VARIABLE DEFINITIONS
%
close all
img_names = {'ref','cor','BiC','CO2','pHe'};
noisenames = {'orig','denoised'};
spnames = {'BiC','CO2'}; %specifies names of spectral bins to define
maps = struct;
roinames = cell(10);
colors = [1 0 0;0 1 0;1 1 0;0 1 1;1 0 1;0 0 1];
pKa = 6.17; %pKa of bicarbonate-CO2
images.CSI.apod = 0; %apodization, in Hz (adjustable using GUI)    
home = pwd;

ovlynamesproc = {['Integral_' spnames{1}],['Integral_' spnames{2}],...
    ['LW_' spnames{1}],['LW_' spnames{2}],['SNR_' spnames{1}],...
    ['SNR_' spnames{2}]}; %used when processing images
for i = 1:length(img_names)-2
    imgname = img_names{i+2};
    ovlynames.(imgname) = {[imgname noisenames{1}],[imgname noisenames{2}]};
end
ovlynames.upsnr = {['UpSNR' spnames{1}],['UpSNR' spnames{2}]};
ovlynames.pctlw = {['pctLW' spnames{1}],['pctLW' spnames{2}]};
ovlynames.diffmaps = {'dpH','rescuedvoxels'};
ovlygrps = fieldnames(ovlynames);

% Define logicals used with plotting (can often toggle with GUI)
%
specflg = false; %used to ID best apodization before denoising
ovlyflg = false; %used to plot overlays on final figure once ROIs defined
gridflg = true; %used to toggle voxel grid on/off
voxflg = false; %used to toggle selected voxel overlay on/off
bndplotflg = false; %used to toggle on/off if integral bounds are displayed
boundflg = false; %if true, user defines different bounds for denoised 
    %dataset (may be necessary if lots of apodization)
alignflg = false; %if true, peaks in CSI will be automatically aligned 
    %across all voxels
finishflg = false; %becomes true when finished with a figure    

diffvox = zeros(2,1); %used for voxel shifting
zoomfac = 0; %factor used in zooming (positive zooms in; negative zooms out)
otp = 0.3; %transparency of 13C overlays 
spax = []; %used for plot axes for spectra

% Check to see if previous images are an input. If so, load images and skip 
% to ROI processing. Otherwise, continue with function as normal.
%
if nargin > 0 && exist('imageprev','var') == 1 && isfield(imageprev,'ref')
    prompt = {'Previous images detected as input. Use?'};
    choices = {'Yes' 'No'};
    answer = listdlg('ListString' , choices , ...
        'SelectionMode' , 'single' , ...
        'ListSize' , [200 30] , ...
        'PromptString' , prompt);
    if answer == 1 
        disp('Loading images from input data...');
        images = imageprev;
        loadflg = false;
        specfig = figure;
        delete(specfig);
    else
        disp('Images will not be loaded from previous data. Continuing...');
        loadflg = true;
    end
else
    loadflg = true;
end


%% DATA PROCESSING: Image Load
%
% Prompt user to specify location of .fid file containing 2D CSI data. 
% Load into structure "images"
%

if loadflg %only perform next 2 sections if images not loaded from previous input
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
%             strings = split(pathname,'/');
%             savename = strings{end-1};
            % Load procpar info from 2D CSI .fid directory
            cd([pathname filename]);
            system([svk_path 'svk_file_convert -i fid -o csixy -t2']);
            ksp = read_ddf_image('csixy');
            ksp.img = permute(ksp.img,[3 2 1]); %make matrix y by x by f
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
            for i = 1:2
                imgname = img_names{i};
                fnames = dir([filestrref{i} '*.fid']);
                filenames = cell(length(fnames),1);
                for j = 1:length(fnames) %read in file names
                    filenames{j} = fnames(j).name;
                end
                if lastflg
                    H1filename = filenames{end};
                else
                    prompt = 'Choose file to load:';
                    choices = filenames;
                    answer = listdlg('ListString' , choices , ...
                             'SelectionMode' , 'single' , ...
                             'ListSize' , [200 50] , ...
                             'PromptString' , prompt);
                    H1filename = filenames{answer};
                end
                disp([imgname ' image: Loading ' H1filename])
                cd(H1filename);
                H1info = load_procpar('procpar');
                %determine which parameters relate to x, y, z; load parameters
                f_id = [{H1info.dimX},{H1info.dimY},{H1info.dimZ}]; %parameter names for xyz-FOV lengths
                o_id = [{H1info.posX},{H1info.posY},{H1info.posZ}]; %parameter names for xyz-offsets
                paramload(f_id,o_id,imgname)
                cd('..')
                [re,im] = load_fids(strtok(H1filename,'.'));
                ksp.(imgname).raw = refimgorder(re,im,images.(imgname).np,...
                    images.(imgname).etl);
                for j = 1:images.(imgname).np(3)
                    images.(imgname).img(:,:,j) = fftshift(fftn(fftshift(...
                        squeeze(ksp.(imgname).raw(:,:,j)))));
                end
                images.(imgname).img = abs(images.(imgname).img);
            end
        end
        cd(home);
        break;
    end
    cd(home);

    % Determine vectors for displaying all images within same FOV + scaling
    % factors
    %
    images.plot.Hxbnds = [-1*images.ref.fov(1),images.ref.fov(1)] ...
        / 2 - images.ref.off(1); %for 1H, axial
    images.plot.Hybnds = [-1*images.ref.fov(2),images.ref.fov(2)] ...
        / 2 - images.ref.off(2);
    images.plot.Kxbnds = [-1*images.cor.fov(1),images.cor.fov(1)] ...
        / 2 - images.cor.off(1); %for 1H, coronal
    images.plot.Kzbnds = [-1*images.cor.fov(2),images.cor.fov(2)] ...
        / 2 - images.cor.off(3);
    images.plot.Pxbnds = [-1*images.CSI.fov(1),images.CSI.fov(1)] ...
        / 2 - images.CSI.off(1); %for pH and associated images
    images.plot.Pybnds = [-1*images.CSI.fov(2),images.CSI.fov(2)] ...
        / 2 - images.CSI.off(2);
    images.plot.Pzbnds = [-1*images.CSI.fov(3),images.CSI.fov(3)] ...
        / 2 - images.CSI.off(3);

    
%% DATA PROCESSING: Zerofilling, Integral/Noise Regions, Grid Shifting/Apodization
%    
    % Prompt user to specify zerofilling for 2D CSI data
    %
    prompt = {'Matrix size, spatial (both directions):' , ...
        'Spectral points:'};
    dlg_title = 'Zero-filling';
    num_lines = 1;
    def = {num2str(images.CSI.npx(1)*2) , num2str(images.CSI.np)};    
    answer2 = inputdlg(prompt , dlg_title , num_lines , def);
    answer2{3} = '0'; %perform spectral deletions, etc w/o apodization
    answer2{4} = '0';
    answer2{5} = '0';
    oldnp = [images.CSI.npx(1),images.CSI.np]; %non-zerofilled dimensions
    newnp = [str2double(answer2{1}),str2double(answer2{2})]; %store 
        %zerofilling values
    answer2{1} = num2str(oldnp(1)); %don't zerofill until denoising
    answer2{2} = num2str(oldnp(2));
    procCSI(answer2,false,false);
    
    images.ppm = linspace(images.CSI.soffp + images.CSI.np/2 * images.CSI.srp,...
        images.CSI.soffp - (images.CSI.np/2 -1) * images.CSI.srp,...
        images.CSI.np); %vector of ppm values to match SIVIC
    bndfig = figure('DeleteFcn',@finishFig); %use for defining integral bounds
    images.spbin = specIntegrals; %get noise region for SNR calculation as
        %well as deletion regions
    finishflg = false; %reset after specifying regions
  
    % Generate images.ppmdelete vector based on spectral deletions, then 
    % convert ppm values for integral bounds to indices for both vectors 
    % (NOTE: ppm will be used with original data, ppmdelete with denoised)
    %
    images.ppmdelete = images.ppm;
    images.ppmdelete(images.spbin.delete) = [];
    if BiCppm(1) < BiCppm(2) %check which is greater, which is lesser
        images.spbin.img.(spnames{1}) = find(images.ppm > BiCppm(1) & ...
            images.ppm < BiCppm(2));
        images.spbin.imgdenoised.(spnames{1}) = find(images.ppmdelete > ...
            BiCppm(1) & images.ppmdelete < BiCppm(2));
    else
        images.spbin.img.(spnames{1}) = find(images.ppm > BiCppm(2) & ...
            images.ppm < BiCppm(1));
        images.spbin.imgdenoised.(spnames{1}) = find(images.ppmdelete > ...
            BiCppm(2) & images.ppmdelete < BiCppm(1));        
    end
    if CO2ppm(1) < CO2ppm(2)
        images.spbin.img.(spnames{2}) = find(images.ppm > CO2ppm(1) & ...
            images.ppm < CO2ppm(2));
        images.spbin.imgdenoised.(spnames{2}) = find(images.ppmdelete > ...
            CO2ppm(1) & images.ppmdelete < CO2ppm(2));        
    else
        images.spbin.img.(spnames{2}) = find(images.ppm > CO2ppm(2) & ...
            images.ppm < CO2ppm(1));
        images.spbin.imgdenoised.(spnames{2}) = find(images.ppmdelete > ...
            CO2ppm(2) & images.ppmdelete < CO2ppm(1));
    end   
    
    procCSI({num2str(oldnp(1)),num2str(oldnp(2)),num2str(images.CSI.origapod), ...
        num2str(diffvox(1)),num2str(diffvox(2))}); %reproc with deletions
    calcMaps; %this will use deletion spectral bins
    
    % Display interactive figure to look at linewidths/lineshapes across
    % image in order to identify best apodization parameters
    %
    nsax = size(images.ref.img,3);
    axsl = round(nsax / 2); 
    specflg = true; 

    ovlyname = ovlynamesproc{1};
end

dispvox = false(images.CSI.npx(1)); %used to display spectra from a single voxel

% Calculate global variables to be used for image display
%
pixratP = images.CSI.npx(1) / (images.CSI.npx(1) + 1);
pixratz = images.CSI.npx(3) / (images.CSI.npx(3) + 2);
offHa = (images.ref.off(3) - images.CSI.off(3)) / ...
    images.ref.ps(3); %difference in slice stack offsets, 1H vs 2D CSI
offHc = (images.CSI.off(2) - images.cor.off(2)) / ...
    images.cor.ps(3); %# of slices the coronal slice stack is offset from pH image
nsHa = round(images.CSI.ps(3) / images.ref.ps(3) / 2) * 2; %number of 1H slices in 13C slice
nsHc = round(images.CSI.fov(2) / images.cor.ps(3) / 2) * 2;
stHa = round((images.ref.np(3) / 2) - (nsHa / 2 - 1) - offHa);%-1;
fnHa = round((images.ref.np(3) / 2) + (nsHa / 2) - offHa);%-1;
stHc = round((images.cor.np(3) / 2) - (nsHc / 2 - 1) - offHc);%-1;
if stHc < 0
    stHc = 0;
end
fnHc = round((images.cor.np(3) / 2) + (nsHc / 2) - offHc);
if fnHc > images.cor.np(3)
    fnHc = images.cor.np(3);
end

if loadflg
    % Display 1H image
    %     
    specfig = figure(100);   
%     figure(66); title(['If function freezes after CONTINUE button pressed, '...
%         'kill this figure']) 
%         %since sometimes the ginput catch requires this figure to be killed    
    plotAxImg;
%     set(specfig,'Position',[0,scrsz(4)*1/6,scrsz(3),scrsz(4)*5/6]*sppi,...
%     'Units','inches');    
    
    % Create uicontrols in figure for selecting voxels, viewing spectra
    %
    bg1 = uibuttongroup('Position',[0 0 1 .15]);
    bg2 = uibuttongroup('Position',[0 .15 .2 .85]);
    %slice selection
    uicontrol(bg1,'Style','text','Position',[0 32 120 20],...
        'String','Active Slice, Axial');
    uicontrol(bg1,'Style','slider','Min',1,'Max',nsax,'Value',axsl,...
        'SliderStep',[1/nsax 3/nsax],'Position',[0 0 100 30],'Callback',...
        @axSelect);
    %select active image data
    uicontrol(bg1,'Style','text','Position',[120 32 120 20],...
        'String','Active Overlay Data');
    uicontrol(bg1,'Style','popupmenu','Position',[120 20 120 10],...
        'String',ovlynamesproc,'Callback',@selOvly);
    %check box for specifying new integral bounds for denoised data, peak
    %alignment
    uicontrol(bg1,'Style','checkbox','Position',[250 30 120 30],'String',...
        'New integral bounds','Callback',@setBounds);
    uicontrol(bg1,'Style','checkbox','Position',[250 0 120 30],'String',...
        'Align peaks','Callback',@setAlign);
    %finish button
    uicontrol(bg1,'Style','pushbutton','Position',[390 0 100 50],...
        'String','CONTINUE','Callback',@finishFig);
    %apodization + voxel shifting button
    uicontrol(bg2,'Style','text','Position',[10 290 80 30],...
        'String','Spectral apodization (Hz):')
    uicontrol(bg2,'Style','edit','Position',[10 260 80 20],...
        'String',num2str(images.CSI.apod),'Callback',@apodData);
    uicontrol(bg2,'Style','pushbutton','Position',[0 220 90 30],...
        'String','Grid shift CSI','Callback',@gridShift);
    %grid/voxel/spectra toggle buttons + slider for voxel selection overlay
%     ib = uicontrol(bg2,'Style','pushbutton','Position',[0 190 90 30],...
%         'String','Turn integral bounds off','Callback',@selIntegralBounds);
    uicontrol(bg2,'Style','text','Position',[0 150 80 40],...
        'String','Overlay Transparency:');
    uicontrol(bg2,'Style','slider','Min',0,'Max',1,'Value',otp,...
        'SliderStep',[0.1 0.3],'Position',[0 130 90 30],'Callback',@setOvlyTransparency);
    uicontrol(bg2,'Style','slider','Min',-0.5,'Max',0.5,'Value',zoomfac,...
        'SliderStep',[0.01 0.1],'Position',[0 90 90 30],'Callback',@setZoom);
    uicontrol(bg2,'Style','pushbutton','Position',[0 60 80 30],...
        'String','Reset zoom','Callback',@resetZoom);
    %update text when running tests, procedures
    ut = uicontrol(bg2,'Style','text','Position',[0 0 80 60],...
        'String','');

    dispVox;
    waitfor(specfig);
    
    
%% DATA PROCESSING: Denoising, pH + Difference Maps Calculation
%
    % Denoise CSI data, then for BOTH original and denoised calculate 
    % integral + SNR + linewidth maps, as well as pHe map from thresholded 
    % BiC and CO2 images and difference maps between original and denoised
    % datasets (NOTE: original dataset will use apodization number 
    % images.CSI.origapod specified at beginning of file!)
    %    
%     images.CSI.imgdenoised = images.CSI.img;
    if alignflg
        disp(['Spectral peaks will be aligned prior to denoising, based on ' ...
            'user specifications.'])
    end 
    procCSI({num2str(newnp(1)) , num2str(newnp(2)) , ...
        num2str(images.CSI.apod) , num2str(diffvox(1)) , ...
        num2str(diffvox(2))},true,alignflg,true); %reprocess + denoise CSI  
        %with deletions, alignment (if specified), and zerofilling
    images.CSI.imgdenoised = images.CSI.img;

    procCSI({num2str(newnp(1)) , num2str(newnp(2)) , ...
        num2str(images.CSI.origapod) , num2str(diffvox(1)) , ...
        num2str(diffvox(2))},false,false,false); %reprocess CSI with specified 
        %apodization and NO spectral deletions and NO spectral alignment 
        %but with zerofilling, for comparison
    disp(['Comparing denoised data with ' num2str(images.CSI.origapod,'%i') ...
        ' Hz-apodized, non-denoised dataset w/o spectral deletions or '...
        'spectral alignment...'])
    % Run calcMaps on denoised data with original spectral bins, just to
    % get SNR map for thresholding. Then check if user indicated he/she  
    % wants to define different integral bounds for denoised data. If not,  
    % keep spectral integral regions for both original and denoised.
    %
    calcMaps(true,true); %calculate for denoised data
    if boundflg 
        bndfig = figure;
        images.spbin.imgdenoised = specIntegrals(true);
    end    
    for i = 1:2
        if i == 1
            calcMaps(false,false); %calculate for non-denoised data
            mapname = 'orig';
            noisename = 'img';
        else
            calcMaps(true,true); %calculate for denoised data
            mapname = 'denoised';
            noisename = 'imgdenoised';
        end
        for j = 1:2 %save BiC and CO2 integral maps as new images
            spname = spnames{j};
            images.([spname mapname]) = maps.(spname).(noisename).sum ...
                .* maps.(spname).(noisename).goodSNR; 
        end
        images.(['pHe' mapname]) = pKa + log10(images.(['BiC' mapname]) ./ ...
            images.(['CO2' mapname]) * images.fcorr);
        images.(['pHe' mapname])(isnan(images.(['pHe' mapname]))) = 0; 
            %replace NaNs with zeroes
        images.(['pHe' mapname])(isinf(images.(['pHe' mapname]))) = 0; 
            %replace Infs with zeroes    
    end
    images.dpH = (images.pHedenoised - images.pHeorig) .* ...
        (images.pHedenoised > 0 & images.pHeorig > 0); %change in pH before and after denoising, 
        %excluding voxels that didn't have good SNR before/after  
    images.rescuedvoxels = (images.pHeorig == 0) & ...
        (images.pHedenoised ~= 0); %map of voxels that have sufficient 
        %SNR as a result of denoising
    for i = 1:2
        spname = spnames{i};
        images.(['UpSNR' spname]) = maps.(spname).imgdenoised.snr ./ ...
            maps.(spname).img.snr .* maps.(spname).imgdenoised.goodSNR; 
            %SNR enhancement factor, so long as SNR is above threshold in 
            %denoised dataset
        LWmean = mean(maps.(spname).img.lw(maps.(spname).img.goodSNR)); 
            %mean linewidth taken over all voxels with SNR > cutoff in
            %ORIGINAL image
        images.(['pctLW' spname]) = (LWmean - maps.(spname).img.lw) ./ ...
            LWmean .* maps.(spname).img.goodSNR * 100;   
            %percent difference from mean of linewidth in original image   
            %across entire SNR-thresholded original image, so long as SNR  
            %was above threshold
        images.maxpctLW.(spname) = max(abs(reshape(images.(['pctLW' spname])...
            ,1,[])));    
    end 
    disp(['Max linewidth % difference from mean: ' num2str(max([...
        images.maxpctLW.(spnames{1}) images.maxpctLW.(spnames{2})]),'%3.0f')]);
end

% Display parameters
%
disp('Image acquisition parameters, 1H:')
disp([num2str(images.ref.np(3)) ' slices prescribed'])
disp([num2str(images.ref.ps(3)) ' mm spacing between slices'])
disp(['axial slice offset = ' num2str(images.ref.off(3)) ' mm'])
disp('Image acquisition parameters, 13C pH:')
disp([num2str(images.CSI.ps(3)) ' mm slice thickness'])
disp(['axial slice offset = ' num2str(images.CSI.off(3)) ' mm'])
disp(['1H slices corresponding to 13C slice: ' num2str(stHa) ...
    ' to ' num2str(fnHa)])
disp(['Bicarb/CO2 images thresholded at SNR > ' num2str(images.pHthr) ...
    ' * noise.'])
disp(['Flip angle correction factor used for calculating pH: ' ...
    num2str(images.fcorr)])


%% DATA PROCESSING: ROI Drawing

% Display all 1H slices in interactive figure. User may specify and edit
% any number of ROIs, which will be displayed on the 1H slices 
% corresponding to the pH slice thickness. 
%
nsax = size(images.ref.img,3);
nscor = size(images.cor.img,3);
axsl = round(nsax / 2);
corsl = round(nscor / 2);
yindP = round(images.CSI.npx(2) / 2 + 0.5 - (images.cor.off(2) ...
    - images.CSI.off(2) + (corsl - images.cor.np(3) / 2 - 0.5) ...
    * images.cor.ps(3)) / images.CSI.ps(2)); 
    %row index of pH data corresponding with active coronal slice
vstp = 0.3; %transparency of voxel selection overlay
endflg = false; %used to finish editing an ROI
voxflg = true; %used to toggle selected voxel overlay on/off
scrsz = get(groot,'ScreenSize');
sppi = get(groot,'ScreenPixelsPerInch');

% Copy over ROI data from previous data, if specified
%
if nargin > 0 %detect if previous ROIs were specified as input to function
    if isfield(imageprev,'AllVoxels')
        disp('Previous ROIs detected as input. Loading data...')
        rois = imageprev;
        nROI = length(fieldnames(imageprev));
        roinames = fieldnames(imageprev);
    elseif isfield(roiprev,'AllVoxels')
        disp('Previous ROIs detected as input. Loading data...')
        rois = roiprev;
        nROI = length(fieldnames(roiprev));        
        roinames = fieldnames(roiprev);
    end
    %remove AllVoxels and AllROIs if present in saved ROIs
    if isfield(rois,'AllVoxels')
        nROI = nROI - 1;
        roinames = roinames(~strcmp(roinames,'AllVoxels'));
    end
    if isfield(rois,'AllROIs')
        nROI = nROI - 1;
        roinames = roinames(~strcmp(roinames,'AllROIs'));
    end
else
    nROI = 0; 
end
actroi = nROI;
nROIctr = nROI; %counts ROIs independent of ROI deletion (so a ROI doesn't accidentally get overwritten)
newname = ['ROI' num2str(nROIctr+1)];

% Display all 1H images
%
disp('Displaying 1H slices for drawing ROIs...')
voxfig = figure('Position',[0,scrsz(4)*1/6,scrsz(3),scrsz(4)*5/6]*sppi,'Units','inches');
a = axes;
c = axes;
curax = a; 
plotAxImg; plotCorImg;

% Create uicontrols in figure for selecting voxels
%
bg1 = uibuttongroup('Position',[0 0 1 .15]);
bg2 = uibuttongroup('Position',[0 .15 .1 .85]);
%slice selection
uicontrol(bg1,'Style','text','Position',[0 30 160 20],...
    'String','Active Slice, Axial');
uicontrol(bg1,'Style','slider','Min',1,'Max',nsax,'Value',axsl,...
    'SliderStep',[1/nsax 3/nsax],'Position',[0 0 140 30],'Callback',@axSelect);
uicontrol(bg1,'Style','text','Position',[140 30 160 20],...
    'String','Active Slice, Coronal');
uicontrol(bg1,'Style','slider','Min',1,'Max',nscor,'Value',corsl,...
    'SliderStep',[1/nscor 3/nscor],'Position',[140 0 140 30],'Callback',@corSelect)
%image selection for ROI drawing
uicontrol(bg1,'Style','text','Position',[280 45 140 15],...
    'String','Image for drawing ROIs:');
uicontrol(bg1,'Style','radiobutton','Position',[300 20 80 20],...
    'String','Axial','Callback',@setAx);
uicontrol(bg1,'Style','radiobutton','Position',[300 0 80 20],...
    'String','Coronal','Callback',@setCor);
%ROI create/edit buttons and inputs
uicontrol(bg1,'Style','text','Position',[440 45 80 15],...
    'String','New ROI name:');
nr = uicontrol(bg1,'Style','edit','Position',[440 20 80 20],...
    'String',newname,'Callback',@nameROI);
uicontrol(bg1,'Style','pushbutton','Position',[540 30 120 30],...
    'String','New ROI','Callback',@newROI);
uicontrol(bg1,'Style','pushbutton','Position',[540 0 120 30],...
    'String','Accept ROI Edits','Callback',@endROIEdits);
%active ROI + delete button
uicontrol(bg1,'Style','text','Position',[700 40 60 15],...
    'String','Active ROI:');
rl = uicontrol(bg1,'Style','popupmenu','Position',[770 45 120 10],...
    'String',roinames,'Callback',@selROI);
uicontrol(bg1,'Style','pushbutton','Position',[900 30 120 30],...
    'String','Edit Active ROI','Callback',@editROI);
uicontrol(bg1,'Style','pushbutton','Position',[700 0 320 30],...
    'String','Delete Active ROI','Callback',@delROI);
%finish button
uicontrol(bg1,'Style','pushbutton','Position',[1060 0 100 60],...
    'String','FINISH ROIS','Callback',@finishROIs);
%grid/voxel toggle buttons + slider for voxel selection overlay
uicontrol(bg2,'Style','text','Position',[0 360 120 40],...
    'String','Voxel Selection Transparency:');
uicontrol(bg2,'Style','slider','Min',0,'Max',1,'Value',vstp,...
    'SliderStep',[0.1 0.3],'Position',[0 330 120 30],'Callback',@setVoxTransparency);
gb = uicontrol(bg2,'Style','pushbutton','Position',[0 280 120 30],...
    'String','Turn grid off','Callback',@selGrid);
vob = uicontrol(bg2,'Style','pushbutton','Position',[0 250 120 30],...
    'String','Turn voxel overlay off','Callback',@selVoxOvly);
uicontrol(bg2,'Style','slider','Min',-0.5,'Max',0.5,'Value',zoomfac,...
    'SliderStep',[0.01 0.1],'Position',[0 200 120 30],'Callback',@setZoom);
uicontrol(bg2,'Style','pushbutton','Position',[0 170 120 30],...
    'String','Reset zoom','Callback',@resetZoom);

waitfor(voxfig);


%% DATA PROCESSING: Statistics Calculations + Display
%
% Extract raw voxel values from each ROI, calculate statistics 
%
for i = 1:nROI %loops through all ROIs to draw
    roiname = roinames{i};
    for j = 3:6 %loop through ovlynames fields
        ovlygrp = ovlygrps{j};
        for k = 1:2 %loop through names within each ovlynames field: pHe, 
            %upsnr, pctlw, diffmaps
            img_name = ovlynames.(ovlygrp){k};
            rois.(roiname).(img_name).img = images.(img_name) .* ...
                rois.(roiname).maskP;
            rois.(roiname).(img_name).raw = ...
                rois.(roiname).(img_name).img(rois.(roiname).(img_name).img ~= 0); 
                %contains all nonzero values
            rois.(roiname).(img_name).raw(isnan(rois.(roiname).(img_name).raw)) ...
                = []; %remove NaNs
            rois.(roiname).(img_name).raw(isinf(rois.(roiname).(img_name).raw)) ...
                = []; %remove Infs
            if ~strcmp(img_name,'rescuedvoxels')
%                 rois.(roiname).(img_name).sum = sum(rois.(roiname).(img_name).raw);
                rois.(roiname).(img_name).mean = mean(rois.(roiname).(img_name).raw);
                rois.(roiname).(img_name).std = std(rois.(roiname).(img_name).raw);
                rois.(roiname).(img_name).max = max(rois.(roiname).(img_name).raw);
                rois.(roiname).(img_name).min = min(rois.(roiname).(img_name).raw);
                disp([roiname ', ' img_name ': Mean = ' ...
                    num2str(rois.(roiname).(img_name).mean,'%10.4f') ' ± ' ...
                    num2str(rois.(roiname).(img_name).std,'%10.4f') ...
                    '. Max = ' num2str(rois.(roiname).(img_name).max,'%10.4f') ...
                    '. Min = ' num2str(rois.(roiname).(img_name).min,'%10.4f') '.']);
%                     '. Sum = ' num2str(rois.(roiname).(img_name).sum,'%10.4f') ...                
            else %display rescued voxels info differently
                rois.(roiname).(img_name).sum = sum(rois.(roiname).(img_name).raw);
                disp([roiname ': ' num2str(rois.(roiname).(img_name).sum,'%i')...
                    ' voxels recovered by denoising']);
            end
        end
    end
end


%% DATA PROCESSING: Image + ROI Display
% Calculate intensity limits for plotting
%
% Generate masks that include selected voxels from all ROIs ('AllROIs') and
% all voxels ('AllVoxels')
%
rois.AllROIs.maskP = false(images.CSI.npx(1));
for i = 1:nROI
    roiname = roinames{i};
    rois.AllROIs.maskP = rois.AllROIs.maskP | rois.(roiname).maskP;
end
rois.AllVoxels.maskP = true(images.CSI.npx(1));

% Set limits for overlay color bars based upon max/min in AllROIs
%
bcmax = max([max(reshape(images.BiCorig .* rois.AllROIs.maskP,1,[])), ...
    max(reshape(images.BiCdenoised .* rois.AllROIs.maskP,1,[])),...
    max(reshape(images.CO2orig .* rois.AllROIs.maskP,1,[])),...
    max(reshape(images.CO2denoised .* rois.AllROIs.maskP,1,[]))]);
if bcmax == 0
    bcmax = 300;
end
% pHlimimg = images.pHe.img .* rois.AllROIs.maskP;
% pHmin = floor(min(min(pHlimimg(pHlimimg ~= 0))) * 10) / 10;% - 0.1;
% pHmax = ceil(max(max(pHlimimg)) * 10) / 10;% + 0.1;
pHmin = 6.7;
pHmax = 8;
if isempty(pHmin) %next two loops should catch an error
    pHmin = 6.5;
end
if isempty(pHmax)
    pHmax = 7.4;
end
bcsnrmax = max([max(reshape(images.(['UpSNR' spnames{1}]) .* ...
    rois.AllROIs.maskP,1,[])) , max(reshape(...
    images.(['UpSNR' spnames{2}]) .* rois.AllROIs.maskP,1,[]))]);
bclwmin = min([min(reshape(images.(['pctLW' spnames{1}]),1,[])) , ...
    min(reshape(images.(['pctLW' spnames{2}]),1,[]))]); %since LW (likely) 
    %matters across entire image, don't base limits on selected voxels in
    %ROIs
bclwmax = max([max(reshape(images.(['pctLW' spnames{1}]),1,[])) , ...
    max(reshape(images.(['pctLW' spnames{2}]),1,[]))]);
dpHmax = max(abs(reshape(images.dpH .* rois.AllROIs.maskP,1,[])));

for i = 1:2
    lims.(img_names{i+2}) = [0 bcmax ; 0 bcmax];
end
lims.pHe = [pHmin pHmax; pHmin pHmax];
lims.upsnr = [0 bcsnrmax; 0 bcsnrmax];
lims.pctlw = [bclwmin bclwmax; bclwmin bclwmax];
lims.diffmaps = [-dpHmax dpHmax; 0 1];

roinames{nROI+1} = 'AllROIs'; 
roinames{nROI+2} = 'AllVoxels';
axsl = round(nsax / 2);
corsl = round(nscor / 2);
vstp = 0.3; 
otp = 0.3; %transparency of 13C overlays
zoomfac = 0;
actroi = 1;
actovlygrp = ovlygrps{1};
ovlyflg = true;
finishflg = false; %may have become true after a figure was closed
bndplotflg = false;

% Create interactive figure which displays selected ROI voxels and all
% overlays at once
% 
finalfig = figure('Position',[0,round(scrsz(4)*1/6),round(scrsz(3)*2/3),round(scrsz(4)*3/4)]...
    *sppi,'Units','inches','DeleteFcn',@finishFig);
plotAxImg; 
plotCorImg; 
plotAxOvly;
bg1 = uibuttongroup('Position',[0 0 1 .1]);
bg2 = uibuttongroup('Position',[0 .1 .1 .9]);
%slice selection
uicontrol(bg1,'Style','text','Position',[0 20 120 20],...
    'String','Active Slice, Axial');
uicontrol(bg1,'Style','slider','Min',1,'Max',nsax,'Value',axsl,...
    'SliderStep',[1/nsax 3/nsax],'Position',[0 10 120 10],'Callback',@axSelect);
uicontrol(bg1,'Style','text','Position',[140 20 120 20],...
    'String','Active Slice, Coronal');
uicontrol(bg1,'Style','slider','Min',1,'Max',nscor,'Value',corsl,...
    'SliderStep',[1/nscor 3/nscor],'Position',[140 10 120 10],'Callback',@corSelect)
%active ROI and overlay
uicontrol(bg1,'Style','text','Position',[300 15 80 15],...
    'String','Active ROI:');
uicontrol(bg1,'Style','popupmenu','Position',[380 20 120 10],...
    'String',roinames,'Callback',@selROI);
uicontrol(bg1,'Style','text','Position',[550 15 80 15],...
    'String','Active overlay group:');
uicontrol(bg1,'Style','popupmenu','Position',[630 20 160 10],...
    'String',fieldnames(ovlynames),'Callback',@selOvlyGrp);
%grid/voxel toggle buttons + slider for voxel selection overlay
uicontrol(bg2,'Style','text','Position',[0 360 120 40],...
    'String','Voxel Selection Transparency:');
uicontrol(bg2,'Style','slider','Min',0,'Max',1,'Value',vstp,...
    'SliderStep',[0.1 0.3],'Position',[0 330 120 30],'Callback',@setVoxTransparency);
uicontrol(bg2,'Style','text','Position',[0 280 120 40],...
    'String','Overlay Transparency:');
uicontrol(bg2,'Style','slider','Min',0,'Max',1,'Value',otp,...
    'SliderStep',[0.1 0.3],'Position',[0 250 120 30],'Callback',@setOvlyTransparency);
gb = uicontrol(bg2,'Style','pushbutton','Position',[0 200 120 30],...
    'String','Turn grid off','Callback',@selGrid);
vob = uicontrol(bg2,'Style','pushbutton','Position',[0 170 120 30],...
    'String','Turn voxel overlay off','Callback',@selVoxOvly);
ib = uicontrol(bg2,'Style','pushbutton','Position',[0 140 120 30],...
    'String','Turn integral bounds off','Callback',@selIntegralBounds);
uicontrol(bg2,'Style','slider','Min',-0.5,'Max',0.5,'Value',zoomfac,...
    'SliderStep',[0.01 0.1],'Position',[0 100 120 30],'Callback',@setZoom);
uicontrol(bg2,'Style','pushbutton','Position',[0 70 120 30],...
    'String','Reset zoom','Callback',@resetZoom);

dispVox;
waitfor(finalfig);


%% DATA PROCESSING: Save Data 
% Write out data to .csv file for easy copying
%
cd(home);
fname = 'mouseproc_output.txt';
fileID = fopen(fname,'w');
fprintf(fileID,'%s \t','ROI name');
fprintf(fileID,'%s \t','pH avg');
fprintf(fileID,'%s \t','pH std');
fprintf(fileID,'%s \t','pH min');
fprintf(fileID,'%s \t','resc vox');
fprintf(fileID,'%s \t','dpH min');
fprintf(fileID,'%s \t','dpH max');
fprintf(fileID,'%s \t','BiC upSNR mean');
fprintf(fileID,'%s \t','BiC upSNR std');
fprintf(fileID,'%s \t','CO2 upSNR mean');
fprintf(fileID,'%s \n','CO2 upSNR std');
for i = 1:nROI
    roiname = roinames{i};
    fprintf(fileID,'%12s \t',roiname);
    data = [rois.(roiname).pHedenoised.mean rois.(roiname).pHedenoised.std ...
        rois.(roiname).pHedenoised.min rois.(roiname).rescuedvoxels.sum ...
        rois.(roiname).dpH.min rois.(roiname).dpH.max rois.(roiname).UpSNRBiC.mean ...
        rois.(roiname).UpSNRBiC.std rois.(roiname).UpSNRCO2.mean ...
        rois.(roiname).UpSNRCO2.std];
    fprintf(fileID,'%10.4f \t %10.4f \t %10.4f\n',data);
end
fclose(fileID);
disp(['ROI statistics written to ' fname ' in directory ' home]);

% Open browser window to save data. Variables saved are images and rois.
%
while 1
    cd(save_dir)
    [sfile , spath] = uiputfile('*.mat' , 'Save Workspace Variables As');
    if isequal(sfile , 0) || isequal(spath , 0)
        prompt = {'Are you sure you do not want to save your data?'};
        choices = {'Yes' 'No'};
        answer4 = listdlg('ListString' , choices , ...
            'SelectionMode' , 'single' , ...
            'ListSize' , [200 30] , ...
            'PromptString' , prompt);
        if answer4 == 1
            cd(home);
            disp('Data were not saved....');
            break;
        end
    else
        save(fullfile(spath , sfile) , 'images' , 'rois');
        cd(home);
        disp(['All images and ROI data saved in ' , fullfile(spath , sfile) , '!']);
        break;
    end
end


%% INTERNAL FUNCTIONS - CALCULATION
%
% calcMaps: Calculate integral, SNR, and linewidth maps. If denoisedflg = 
% true, denoised rather than original data are used. If deleteflg = true,
% spectral bins with deletions are used. If shiftcalcflg = true, update
% shiftmap used for aligning spectra
%
function calcMaps(denoisedflg,deleteflg,shiftcalcflg)
if nargin < 1
    deleteflg = true;
    denoisedflg = false;
    shiftcalcflg = false;
elseif nargin < 2
    deleteflg = true;    
    shiftcalcflg = false;
elseif nargin < 3
    shiftcalcflg = false;
end
if denoisedflg
    nname = 'imgdenoised'; %use denoised data
else
    nname = 'img'; %use non-denoised data
end
if deleteflg
    binname = 'imgdenoised';
else
    binname = 'img';
end
% % Calculate average noise standard deviation in all voxels. To correct for 
% % any linear baseline slope, split noise region in half, reverse latter
% % half, and add to former half. Then take std, multiply by sqrt(2) to
% % account for noise reduction by summation
% halfpt = floor(length(images.spbin.noise.(binname))/2);
% sumnoise = images.CSI.(nname)(:,:,images.spbin.noise.(binname)(1:halfpt)) + ...
%     flip(images.CSI.(nname)(:,:,images.spbin.noise.(binname)(halfpt+1:2*halfpt))...
%     ,3);
% noisemap = abs(std(sumnoise/2,0,3)) * sqrt(2);   
% Calculate average noise standard deviation in all voxels
noisemap = abs(std(images.CSI.(nname)(:,:,images.spbin.noise.(binname)),0,3));   
% noise = reshape(images.CSI.(nname)(:,:,images.spbin.noise),1,[]);
% images.CSI.noise.(imgname) = abs(std(noise)); %takes from all voxels
% Calculate maps of integral, SNR, linewidth over all images
for jj = 1:length(spnames)
    spn = spnames{jj};
    maps.(spn).(nname).img = abs(...
        images.CSI.(nname)(:,:,images.spbin.(binname).(spn))); 
        %select spectral range to process, take abs val 
    maps.(spn).(nname).sum = sum(maps.(spn).(nname).img,3);
    maps.(spn).(nname).max = max(maps.(spn).(nname).img,[],3);
    maps.(spn).(nname).snr = maps.(spn).(nname).max ./ noisemap;
    %Calculate linewidths
    hm = maps.(spn).(nname).max / 2; %half-max value for each voxel
    for kk = 1:images.CSI.npx(1)
        for ll = 1:images.CSI.npx(2)
            if maps.(spn).(nname).snr(ll,kk) > images.pHthr*2/3 %only calculate
                %LW if SNR is high enough; use 1/2 the threshold for shift
                %analysis
                maxind = find(abs(images.CSI.(nname)(ll,kk,:)) == ...
                    maps.(spn).(nname).max(ll,kk)); %index of voxel max value
                maxind = maxind(maxind >= images.spbin.(binname).(spn)(1) & ...
                    maxind <= images.spbin.(binname).(spn)(end)); %make sure 
                    %>1 value wasn't caught; make sure it's in spectral bin
                hmind = zeros(2,1);
                for mm = 1:round(3*length(images.spbin.(binname).(spn)))   
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
                maps.(spn).(nname).lw(ll,kk) = (hmind(2) - hmind(1)) ...
                    * images.CSI.sr; %peak linewidth, in Hz
                maps.(spn).(nname).midpt(ll,kk) = mean(hmind); %midpoint 
                    %index inbetween half-max values (used for shift
                    %correction)
            else
                maps.(spn).(nname).lw(ll,kk) = 0;
                maps.(spn).(nname).midpt(ll,kk) = 0;
            end
        end
    end
    maps.(spn).(nname).goodSNR = maps.(spn).(nname).snr > ...
        images.pHthr; %used to display only good-SNR voxels
    %Save images as overlays accessible to plotting
    mname = ['Integral_' spn];
    maps.(mname) = maps.(spn).(nname).sum; 
    lims.(mname) = [min(reshape(maps.(mname),1,[])) ...
        max(reshape(maps.(mname),1,[]))];
    mname = ['SNR_' spn];
    maps.(mname) = maps.(spn).(nname).snr; 
    lims.(mname) = [0 max(reshape(maps.(mname),1,[]))];
    mname = ['LW_' spn];
    maps.(mname) = maps.(spn).(nname).lw; 
    lims.(mname) = [0 max(reshape(maps.(mname) .* ...
        maps.(spn).(nname).goodSNR,1,[]))];            
end
%Calculate map of spectral shifts required (in indices rather than ppm)
%in order to align peaks, as follows. First, calculate midpoint 
%inbetween peaks for each spectrum. Then take mean value and use to 
%generate map of shifts for each spectra (to be used as input to
%procCSI)
if shiftcalcflg %calculate map for shifting spectra to align (used in procCSI)
%     midptmap = (maps.BiC.(nname).midpt + maps.CO2.(nname).midpt) ./ 2;
    midptmap = maps.BiC.(nname).midpt; %use BiC only
%     mpmean = mean(midptmap((maps.BiC.(nname).snr >= (images.pHthr*2/3)) & ...
%         (maps.CO2.(nname).snr >= (images.pHthr*2/3)))); %take mean of all midpoints 
%         %where both metabolites have SNR above 2/3 the cutoff
    mpmean = mean(midptmap((maps.BiC.(nname).snr >= (images.pHthr*2/3)))); 
        %take mean of all midpoints where BiC has SNR above 2/3 the cutoff        
    images.shiftmap = midptmap - mpmean; %map of shift values
%     images.shiftmap((maps.BiC.(nname).snr < (images.pHthr*2/3)) | ...
%         (maps.CO2.(nname).snr < (images.pHthr*2/3))) = 0; %don't shift spectra 
%         %with bad SNR, either metabolite
    images.shiftmap((maps.BiC.(nname).snr < (images.pHthr*2/3))) = 0;  
        %don't shift spectra with bad BiC SNR        
    images.shiftmap(abs(images.shiftmap) > 2.5) = 0; %don't shift if it's 
        %above threshold value    
end
%Adjust plotting limits so both metabolite maps have same scaling
lims.(['Integral_' spnames{1}])(2) = max(...
    [lims.(['Integral_' spnames{1}])(2) lims.(['Integral_' spnames{2}])(2)]);
lims.(['Integral_' spnames{2}]) = lims.(['Integral_' spnames{1}]);
lims.(['SNR_' spnames{1}])(2) = max(...
    [lims.(['SNR_' spnames{1}])(2) lims.(['SNR_' spnames{2}])(2)]);
lims.(['SNR_' spnames{2}]) = lims.(['SNR_' spnames{1}]);
lims.(['LW_' spnames{1}])(2) = max(...
    [lims.(['LW_' spnames{1}])(2) lims.(['LW_' spnames{2}])(2)]);
lims.(['LW_' spnames{2}]) = lims.(['LW_' spnames{1}]);
end


%% INTERNAL FUNCTIONS - IMAGE/SPECTRAL PROCESSING
%
% paramload: Takes loaded procpar info, saves into images._
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

% specIntegrals: Prompts user to specify spectral regions for calculation,
% either for peaks or for noise + deletion regions depending on input 
% logical spbinsflg
%
function bins = specIntegrals(spbinsflg)
% Plot all spectra (from newly apodized image), prompt user to 
% specify spectral bins of interest and/or noise 
%
if nargin < 1 %by default, only calculate noise
    spbinsflg = false;
end
figure(bndfig); 
clf(bndfig); 
hold on
for ii = 1:images.CSI.npx(1)
    for jj = 1:images.CSI.npx(2)
        if isfield(maps,'BiC') %only plot voxels with good SNR
            if maps.BiC.imgdenoised.goodSNR(ii,jj) && ...
                maps.CO2.imgdenoised.goodSNR(ii,jj)
                nspec = squeeze(abs(images.CSI.imgdenoised(ii,jj,:)) ...
                    / max(reshape(abs(images.CSI.imgdenoised),1,[])));
                plot(images.ppmdelete,nspec); set(gca,'xdir','reverse');                
            end
        else
            nspec = squeeze(abs(images.CSI.img(ii,jj,:)) ...
                / max(reshape(abs(images.CSI.img),1,[])));
            plot(images.ppm,nspec); set(gca,'xdir','reverse');
        end
    end
end
pause(.1); %otherwise, won't plot spectra for some reason...
figure(bndfig);
if spbinsflg %either get spectral bins, or get noise bin; not both!
    for ii = 1:length(spnames)
        spname = spnames{ii};
        title([spname ': Select spectral region'])
        nb = getrect;
        bins.(spname) = find(images.ppmdelete > nb(1) & ...
            images.ppmdelete < (nb(1) + nb(3)));
        yvec = [0 1];
        xvec = [1 1] * images.ppmdelete(bins.(spname)(1));
        line(xvec,yvec,'Color',colors(ii,:),'LineStyle','--')
        xvec = [1 1] * images.ppmdelete(bins.(spname)(end));
        line(xvec,yvec,'Color',colors(ii,:),'LineStyle','--')
    end
else
    title('Select noise region')
    nb = getrect;
    bins.noise.img = find(images.ppm > nb(1) & images.ppm < (nb(1) + nb(3)));
    yvec = [0 1];
    xvec = [1 1] * images.ppm(bins.noise.img(1));
    line(xvec,yvec,'Color',colors(3,:),'LineStyle','--')
    xvec = [1 1] * images.ppm(bins.noise.img(end));
    line(xvec,yvec,'Color',colors(3,:),'LineStyle','--')
    title('Select deletion regions. When finished, close figure')
    bins.delete = [];
    while ~finishflg
        try
            nb = getrect;
            yvec = [1 ; 1];
            xvec = [nb(1) (nb(1) + nb(3))];
            area(xvec,yvec,'FaceAlpha',0.2,'EdgeAlpha',0.4,...
                'LineWidth',0.1,'FaceColor','black');
            bins.delete = [bins.delete find(images.ppm > nb(1) & ...
                images.ppm < (nb(1) + nb(3)))];
        catch
        end
    end
    %When done with deletion regions, shift over noise indices based on 
    %what was deleted
    bins.noise.imgdenoised = bins.noise.img - sum(bins.delete < ...
        bins.noise.img(1));
end
end

% procCSI: Uses SIVIC scripts to reprocess and denoise raw CSI data. If 
% deleteflg = true, remove spectral regions at end of processing. If 
% shiftflg = true, shift over each spectrum in order to align BiC + CO2 
% peaks with one another (shifting involves both k-space shifting and 
% moving points over). If denoiseflg = true, denoise CSI data using triple
% Tucker method.
%
function procCSI(parsin,deleteflg,shiftflg,denoiseflg)
if nargin < 2
    shiftflg = false;
    deleteflg = true;
    denoiseflg = false;
elseif nargin < 3
    shiftflg = false;
    denoiseflg = false;
elseif nargin < 4
    denoiseflg = false;
end
% cd([pathname filename]);
% system([svk_path 'svk_apodize -i fid -o fidx -t4 -f1 --width ' ...
%      parsin{3}]);
% system([svk_path 'svk_zerofill -i fidx.dcm -o fidx -t4 --custom ' parsin{2}]);
%     %this zerofills spectral only;    
% system([svk_path 'svk_zerofill -i fidx.dcm -o fidx -t4 --custom ' ...
%      parsin{2} ' --customx ' parsin{1} ' --customy ' parsin{1}]);
% % system([shift_path ' ' -answer{4} ' ' answer{5}]); %x-shift needs to be 
%     %flipped here, probably due to SIVIC flipping
% system([svk_path 'svk_fft -i fidx.dcm -o csixy -t2 --vsx ' parsin{4} ...
%     ' --vsy ' parsin{5}]); %x-shift needs to be flipped here, probably due 
%     %to SIVIC flipping   
% system([svk_path 'svk_fft -i fid -o csixy -t2 --vsx ' parsin{4} ...
%     ' --vsy ' parsin{5}]); %x-shift needs to be flipped here, probably due 
%     %to SIVIC flipping 
% system([svk_path 'svk_file_convert -i fid -o csixy -t2']);
% temp = read_ddf_image('csixy');
% cd(home);
temp = ksp.img;
temp = flip(temp,1);
% remove deletion regions specified in spbin field (if deleteflg)
temp = fftshift(fftn(temp)); %to spectrum
if deleteflg 
    temp(:,:,images.spbin.delete) = [];
end
% shifting over points (if shiftflg = true)
if shiftflg
    ishft = round(images.shiftmap); %map of integer shifts (in indices) 
        %that will be performed via circshift
    kshft = images.shiftmap - ishft; %map of shifts (in indices) that will
        %be performed using FID multiplication, ranging from +/-0.5
    kshft = kshft * ksp.ddf.sweepwidth / ksp.ddf.specpoints; 
        %convert kshft to units of frequency       
    for ii = 1:size(temp,1)
        for jj = 1:size(temp,2)
            temp(jj,ii,:) = circshift(temp(jj,ii,:),...
                -ishft(jj,ii)); %ishft needs to be negative for best results
        end
    end
end
temp = ifftn(temp); %back to FID
dt = ksp.ddf.dwelltime/1000;
t(1,1,:) = 0:dt:dt*(size(temp,3)-1); %time vector for apodization/
    %shifting 
%spectral shifting via complex exponential (if shiftflg = true)
if shiftflg 
    shftmat = zeros(size(temp));
    for ii = 1:size(temp,1)
        for jj = 1:size(temp,2)
            shftmat(jj,ii,:) = exp(-1i * 2 * pi * t * kshft(jj,ii));
        end
    end
    temp = fft(fft(temp,[],1),[],2); %FT in spatial 
        %dimensions so shift factors line up w/ each spatial voxel
    temp = temp .* shftmat; 
    temp = ifft(ifft(temp,[],1),[],2); %inverse FT  
        %in spatial domain
end
%spatial grid shifting - MAY NEED TESTING!
kmax = 1 / 2 / ksp.ddf.pixel_spacing(1); %k-max, in 1/mm
kvals = linspace(-kmax,kmax,ksp.ddf.npix(1));
shiftx = (str2double(parsin{4})-.5) * ksp.ddf.pixel_spacing(1) * ...
    ksp.ddf.npix(1) / str2double(parsin{1}); %spatial x-shift, in mm
shifty = (str2double(parsin{5})-.5) * ksp.ddf.pixel_spacing(2) * ...
    ksp.ddf.npix(2) / str2double(parsin{1}); %spatial y-shift, in mm
imshft = exp(1i * 2 * pi * kvals' * shifty) * ...
    exp(1i * 2 * pi * kvals * shiftx);
imshft = repmat(imshft, [1 1 size(temp,3)]);
temp = temp .* imshft;
%apodization (+ progressive denoising if denoiseflg = true)
if denoiseflg
    nDN = ceil(images.CSI.apod / images.CSI.apodincr); %number of  
        %apodization-denoising cycles
    progapod = str2double(parsin{3}) / nDN; %apodization used per cycle, in Hz
    sigma = 1 * sqrt(2 * log(2)) / pi / progapod; %sigma parameter of 
        %Gaussian
    apdmat = repmat(exp(-1 * (t .^ 2 / 2 / (sigma ^ 2))),ksp.ddf.npix(1), ...
        ksp.ddf.npix(2),1); %gaussian            
    for ii = 1:nDN
        temp = temp .* apdmat; %apodize
        temp = tensorDenoising(temp,sizeLR); %denoise
    end
    % deapdmat = repmat(exp(t2*images.CSI.apod*pi),ksp.ddf.npix(1), ...
    %     ksp.ddf.npix(2),1); %exponential; an additional factor of 2*pi is 
    %     %required to match SIVIC apodization results
    deapdmat = repmat(exp(nDN * (t .^ 2 / 2 / (sigma ^ 2))),ksp.ddf.npix(1), ...
        ksp.ddf.npix(2),1); %gaussian     
    temp = temp .* deapdmat;
    temp = fftn(temp); %back to spectra
    temp = tensorDenoising(temp,sizeLR); %final denoise
    % Re-apodize data to match apodization of original dataset
    temp = ifftn(temp); %to FID
    % reapdmat = repmat(exp(-t2*images.CSI.origapod*pi),ksp.ddf.npix(1), ...
    %     ksp.ddf.npix(2),1); %exponential
    sigma = 1 * sqrt(2 * log(2)) / pi / images.CSI.origapod; %sigma parameter of 
        %Gaussian
    reapdmat = repmat(exp(-1 * (t .^ 2 / 2 / (sigma ^ 2))),ksp.ddf.npix(1), ...
        ksp.ddf.npix(2),1); %gaussian 
    temp = temp .* reapdmat;
else
%     apdmat = repmat(exp(-t*str2double(parsin{3})*pi),ksp.ddf.npix(1), ...
%         ksp.ddf.npix(2),1); %exponential
    sigma = 1 * sqrt(2 * log(2)) / pi / str2double(parsin{3}); %sigma parameter 
        %of Gaussian
    apdmat = repmat(exp(-1 * (t .^ 2 / 2 / (sigma ^ 2))),ksp.ddf.npix(1), ...
        ksp.ddf.npix(2),1); %gaussian 
    temp = temp .* apdmat; %apodize    
end
%zerofilling
pad = (str2double(parsin{1}) - ksp.ddf.npix(1)) / 2;
temp = padarray(temp,[1 1 0]*pad); %zerofill spatially
pad = (str2double(parsin{2}) - ksp.ddf.specpoints);
temp = padarray(temp,[0 0 1]*pad,0,'post'); %zerofill spectrally
images.CSI.img = fftn(temp); %to spectra
% images.CSI.img = flip(images.CSI.img,1); 
% Update CSI parameters based upon zerofilling + loading in image
images.CSI.fov = [info.lpe info.lpe2 (info.thk / 10)] * 10;
    %FOV, in mm: 1st PE is x, 2nd PE is y
images.CSI.sw = info.sw; %spectral width, in Hz
images.CSI.npx = [str2double(parsin{1}) str2double(parsin{1}) 1]; 
    %# of pixels
images.CSI.ps = images.CSI.fov ./ images.CSI.npx;
    %#pixel size, in mm
images.CSI.off = [info.ppe+globshft(1)-str2double(parsin{4})*images.CSI.ps(1) ...
    info.ppe2+globshft(2)-str2double(parsin{5})*images.CSI.ps(2) info.pss0];  
images.CSI.soffp = ksp.ddf.ppm_ref; %spectral center, in ppm
images.CSI.np = str2double(parsin{2}); %# of spectral points
% if isfield(images,'spbin') %adjust number of spectral points due to deletion
%     if isfield(images.spbin,'delete')
%         images.CSI.np = images.CSI.np - length(images.spbin.delete);
%     end
% end
images.CSI.sr = images.CSI.sw / images.CSI.np; %spectral resolution, in Hz
images.CSI.srp = images.CSI.sr / ksp.ddf.transmit_freq; %spectral 
    %resolution, in ppm
images.CSI.apod = str2double(parsin{3}); %spectral apodization, in Hz
%Adjust some of the plotting parameters
images.plot.Pxbnds = [-1*images.CSI.fov(1),images.CSI.fov(1)] ...
    / 2 + images.CSI.off(1); %for 13C
images.plot.Pybnds = [-1*images.CSI.fov(2),images.CSI.fov(2)] ...
    / 2 + images.CSI.off(2);
images.plot.Pzbnds = [-1*images.CSI.fov(3),images.CSI.fov(3)] ...
    / 2 + images.CSI.off(3);
end

% apodData: Adjusts apodization, apodizes 2D CSI data in spectral dimension
% and recalculates all images
%
function apodData(source,~)
images.CSI.apod = str2double(source.String);
if alignflg %need to proc CSI w/o shifts, calculate shiftmap, then apply to 
    %data
    procCSI({num2str(oldnp(1)),num2str(oldnp(2)),num2str(images.CSI.apod), ...
        num2str(diffvox(1)),num2str(diffvox(2))},true,false);
    calcMaps(false,true,true);
end
procCSI({num2str(oldnp(1)),num2str(oldnp(2)),num2str(images.CSI.apod), ...
    num2str(diffvox(1)),num2str(diffvox(2))},true,alignflg);
calcMaps; 
plotAxImg;
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
diff = -[(xpos - xposvox),(ypos - yposvox)];
diffvox = diff ./ images.CSI.ps(1:2); %convert to # of voxels
%Reprocess CSI data with voxel shifts, plot
if alignflg %need to proc CSI w/o shifts, calculate shiftmap, then apply to 
    %data
    procCSI({num2str(oldnp(1)),num2str(oldnp(2)),num2str(images.CSI.apod), ...
        num2str(diffvox(1)),num2str(diffvox(2))},true,false);
    calcMaps(false,true,true);
end
procCSI({num2str(oldnp(1)) , num2str(oldnp(2)) , ...
    num2str(images.CSI.apod) , num2str(diffvox(1)) , num2str(diffvox(2))},...
    true,alignflg);
calcMaps; 
plotAxImg;
end


%% INTERNAL FUNCTIONS - IMAGE/SPECTRUM DISPLAY
%
% plotSpecVox: Plots CSI spectra from selected voxel, along with integral
% region bounds (depending on logical bndplotflg). Input logical origflg
% indicates if original or denoised bounds should be used.
%
function plotSpecVox(bndplotflg)
if nargin < 1
    bndplotflg = false;
end
delete(spax); %clear currently plotted spectra
clear spax
if ~isempty(dispvox(dispvox)) %check to see if any voxels selected for plotting
    [sprow,spcol] = find(dispvox); %get coordinates for spectra to plot
    if ovlyflg
        spxl = max(spcol) - min(spcol) + 1;
        spyl = max(sprow) - min(sprow) + 1;
        xl = .4 / max([spxl/2 spyl]); %horizontal length of one spectrum plot
        yl = xl; %vertical length of one spectrum plot
        xp = .5 + xl * (spcol - min(spcol)) / 2;
        yp = .95 - yl * (sprow - min(sprow) + 1);
        yp((length(sprow) + 1):(length(sprow)*2)) = ...
            .5 - yl * (sprow - min(sprow) + 1);
    else
        xl = .25;
        yl = .35;
        xp = .6 + .5 * xl;
        yp = .9 - 1.5 * yl;
    end
    if ovlyflg
        annotation('textbox',[.7 .9 .3 .1],'String','Original',...
            'EdgeColor','none','FontSize',16)
        annotation('textbox',[.7 .435 .3 .1],'String','Denoised',...
            'EdgeColor','none','FontSize',16)        
        lims.plot = [0 max(abs(reshape(images.CSI.img(sprow,spcol,:)...
            ,1,[]))) * 1.1];
    else
        lims.plot = [0 max(abs(reshape(images.CSI.img,1,[]))) * 1.1];
    end
    for jj = 1:length(sprow)
        %plot original spectrum first
        spax(jj) = axes('Position',[xp(jj) yp(jj) xl yl]);
        plot(abs(squeeze(images.CSI.img(sprow(jj),spcol(jj),:))));
        axis('square'); 
        axis([1 images.CSI.np lims.plot]);
        set(gca,'Xtick',[]); 
        set(gca,'Ytick',[]);
        if bndplotflg %plot integral bounds if specified
            hold on
            for ii = 1:length(spnames)
                spname = spnames{ii};
                imgname = 'img';
                yvec = [1 ; 1] * lims.plot(2);
                xvec = [images.spbin.(imgname).(spname)(1) ...
                    images.spbin.(imgname).(spname)(end)];
                area(xvec,yvec,'FaceAlpha',0.2,'EdgeAlpha',0.4,...
                    'LineWidth',0.1,'FaceColor',colors(ii,:));
            end
            hold off
        end
        %then plot denoised below
        if ovlyflg
            spax(jj+length(sprow)) = axes('Position',...
                [xp(jj) yp(jj+length(sprow)) xl yl]);
            plot(abs(squeeze(images.CSI.imgdenoised(sprow(jj),spcol(jj),:))));
            axis('square'); 
            lims.plot = [0 max(abs(reshape(...
                images.CSI.imgdenoised(sprow,spcol,:),1,[]))) * 1.1];
            axis([1 size(images.CSI.imgdenoised,3) lims.plot]);
            set(gca,'Xtick',[]); 
            set(gca,'Ytick',[]);
            lims.plot = [0 max(abs(reshape(...
                images.CSI.img(sprow,spcol,:),1,[]))) * 1.1]; %return to 
                %original plotbounds for next spectral box     
            if bndplotflg %plot integral bounds if specified
                hold on
                for ii = 1:length(spnames)
                    spname = spnames{ii};
                    imgname = 'imgdenoised';
                    yvec = [1 ; 1] * lims.plot(2);
                    xvec = [images.spbin.(imgname).(spname)(1) ...
                        images.spbin.(imgname).(spname)(end)];
                    area(xvec,yvec,'FaceAlpha',0.2,'EdgeAlpha',0.4,...
                        'LineWidth',0.1,'FaceColor',colors(ii,:));                                       
                end
                hold off
            end    
        end
    end
else
    spax = [];
end
end

% dispVox: Allows user to select voxel for seeing spectra
%
function dispVox(~,~)
while 1
    drawnow; %update variables
    if ~finishflg 
        if isvalid(specfig) %this loop breaks when finishflg is 
        %true or if GUI figure has been deleted    
            try
                figure(specfig);
                [xpos,ypos] = ginput(1); %select 1 voxel
            catch
            end
        elseif ovlyflg
            try
                figure(finalfig);
                nb = getrect(a);
                xx = nb(1):images.CSI.ps(1):sum([nb(1),nb(3)]);
                yy = nb(2):images.CSI.ps(1):sum([nb(2),nb(4)]);
                xpos = sort(repmat(xx',length(yy),1));
                ypos = repmat(yy',length(xx),1);
            catch 
            end
        end 
    end
    drawnow; %update variables
    if finishflg %this loop breaks when finishflg is 
        %true or if GUI figure has been deleted
        return;
    elseif exist('xpos','var')
        % ID selected pH voxels on figures (only for axial currently)
        dispvox = false(images.CSI.npx(1));
        for m = 1:images.CSI.npx(1)
            xl = images.plot.Pxbnds(1) + images.CSI.ps(1) * (m - 1);
            xu = images.plot.Pxbnds(1) + images.CSI.ps(1) * m;
            for n = 1:images.CSI.npx(2)
                yl = images.plot.Pybnds(1) + images.CSI.ps(1) * (n - 1);
                yu = images.plot.Pybnds(1) + images.CSI.ps(1) * n;
                for o = 1:length(xpos)
                    if xpos(o) > xl && xpos(o) < xu && ypos(o) > yl && ...
                        ypos(o) < yu 
                        dispvox(n,m) = true;
                        xpos(o) = []; %remove entry to speed up future runs
                        ypos(o) = [];
                        break;
                    end
                end
            end
        end
    % Plot selected CSI voxel(s) on the other half of the figure
    plotSpecVox(bndplotflg);
    if ovlyflg %update images with outline boxes around selected spectra
        plotAxImg;
        plotAxOvly;
    end
    end
end
end

% plotAxImg: plots axial 1H image (+ grid corresponding to pH voxels, if
% within 13C pH slab)
%
function plotAxImg()
    % First, convert image to grayscale (necessary for overlays)
    plotimg = squeeze(images.ref.img(:,:,axsl));
    hlims = [min(min(plotimg)) max(max(plotimg))];
    plotimg = repmat(mat2gray(plotimg,hlims),[1 1 3]);
    if specflg
        figure(specfig); 
        subplot(1,7,2:5);
    elseif ~ovlyflg
        a = subplot(1,2,1); 
    else
        a = subplot(2,4,1);
    end
    imagesc(images.plot.Hxbnds,images.plot.Hybnds,squeeze(plotimg)); 
    title(num2str(axsl)); 
    axis('equal','off');
    hold on
    if axsl >= stHa && axsl <= fnHa %check if current slice is within 13C pH slab
        % First plot grid
        if gridflg
            for l = images.plot.Pxbnds(1):images.CSI.ps(1):images.plot.Pxbnds(2)
                x = [l l];
                y = [images.plot.Pybnds(1) images.plot.Pybnds(2)];
                plot(x,y,'-g')
            end
            for l = images.plot.Pybnds(1):images.CSI.ps(2):images.plot.Pybnds(2)
                x = [images.plot.Pxbnds(1) images.plot.Pxbnds(2)]; 
                y = [l l];  
                plot(x,y,'-g')
            end
        end
        % Then plot either spectra-related overlay (specfig) or ROI voxels
        % (voxfig)
        if specflg
            v = imagesc((images.plot.Pxbnds-images.CSI.off(1))*pixratP+images.CSI.off(1),...
                (images.plot.Pybnds-images.CSI.off(2))*pixratP+images.CSI.off(2),...
                maps.(ovlyname),lims.(ovlyname));
            colormap(jet); 
            colorbar; 
            axis('square');
            if ~isempty(strfind(ovlyname,spnames{1})) %show only voxels with
                %good SNR for the specific metabolite 
                set(v,'AlphaData',maps.(spnames{1}).img.goodSNR*otp);
            elseif ~isempty(strfind(ovlyname,spnames{2}))
                set(v,'AlphaData',maps.(spnames{2}).img.goodSNR*otp);
            end
        elseif voxflg
            for l = 1:nROI
                name = roinames{l};
                cimg = cat(3,ones(images.CSI.npx(1))*colors(l,1),...
                    ones(images.CSI.npx(1))*colors(l,2),...
                    ones(images.CSI.npx(1))*colors(l,3));
                v = imagesc(images.plot.Pxbnds*pixratP,images.plot.Pybnds*pixratP,...
                    cimg);
                set(v,'AlphaData',rois.(name).maskP*vstp);
            end
            if ovlyflg && sum(sum(dispvox)) > 0 %outline selected voxels
                [olvoxy,olvoxx] = find(dispvox);
                xl = images.plot.Pxbnds(1) + images.CSI.ps(1) * (min(olvoxx)-1);
                xu = images.plot.Pxbnds(1) + images.CSI.ps(1) * max(olvoxx);
                yl = images.plot.Pybnds(1) + images.CSI.ps(1) * (min(olvoxy) - 1);
                yu = images.plot.Pybnds(1) + images.CSI.ps(1) * max(olvoxy);
                plot([xl xu xu xl xl],[yl yl yu yu yl],'g-','LineWidth',2);
            end
        end
    end
    % Plot and label scale bar   
    plotaxes = [images.plot.Hxbnds,images.plot.Hybnds] .* ones(1,4) ...
        * (1-zoomfac);
    barx = [plotaxes(2)*0.9-4 plotaxes(2)*0.9];
    bary = [plotaxes(4) plotaxes(4)] * 0.9;
    plot(barx,bary,'w-','LineWidth',4)
    text(mean(barx),bary(1)*.8/.9,'4 mm','FontSize',12,'Color','w',...
        'HorizontalAlignment','center');
    axis(plotaxes); %adjust zoom
    hold off
end

% plotCorImg: plots coronal 1H image (+ grid corresponding to pH voxels, if
% within 13C pH slab)
%
function plotCorImg()
    % First, convert image to grayscale (necessary for overlays)
    plotimg = squeeze(images.cor.img(:,:,corsl));
    hlims = [min(min(plotimg)) max(max(plotimg))];
    plotimg = repmat(mat2gray(plotimg,hlims),[1 1 3]);
    if ~ovlyflg
        c = subplot(1,2,2);
    else
        c = subplot(2,4,2);
    end
    imagesc(images.plot.Kxbnds,images.plot.Kzbnds,squeeze(plotimg)); 
    title(num2str(corsl)); axis('equal','off')
    hold on
    if corsl >= stHc && corsl <= fnHc %check if current slice is within 13C pH FOV
        % First plot grid
        if gridflg
            for l = images.plot.Pxbnds(1):images.CSI.ps(1):images.plot.Pxbnds(2)
                x = [l l];
                y = [images.plot.Pzbnds(1) images.plot.Pzbnds(2)];
                plot(x,y,'-g')
            end
            for l = images.plot.Pzbnds(1):images.CSI.ps(3):images.plot.Pzbnds(2)
                x = [images.plot.Pxbnds(1) images.plot.Pxbnds(2)];
                y = [l l];  
                plot(x,y,'-g')
            end
        end
        % Then plot ROI voxels
        if voxflg
            for l = 1:nROI
                name = roinames{l};
                cimg = cat(3,ones(1,images.CSI.npx(1))*colors(l,1),...
                    ones(1,images.CSI.npx(1))*colors(l,2),...
                    ones(1,images.CSI.npx(1))*colors(l,3));
                v = imagesc(images.plot.Pxbnds*pixratP,images.plot.Pzbnds*pixratz,...
                    cimg);
                set(v,'AlphaData',squeeze(rois.(name).maskP(yindP,:))*vstp);
            end
        end
    end
    % Plot and label scale bar   
    plotaxes = [images.plot.Kxbnds,images.plot.Kzbnds] .* ones(1,4) ...
        * (1-zoomfac);
    barx = [plotaxes(2)*0.9-4 plotaxes(2)*0.9];
    bary = [plotaxes(4) plotaxes(4)] * 0.9;
    plot(barx,bary,'w-','LineWidth',4)
    text(mean(barx),bary(1)*.8/.9,'4 mm','FontSize',12,'Color','w',...
        'HorizontalAlignment','center');
    axis(plotaxes); %adjust zoom
    hold off
end

% plotAxOvly: plots axial 1H image (+ overlay voxels, if within 13C pH 
% slab)
%
function plotAxOvly()
    % First, convert image to grayscale (necessary for overlays)
    plotimg = squeeze(images.ref.img(:,:,axsl));
    hlims = [min(min(plotimg)) max(max(plotimg))];
    plotimg = repmat(mat2gray(plotimg,hlims),[1 1 3]);
    for l = 1:2 %cycle through all overlays specified in active ovlynames group 
        oname = ovlynames.(actovlygrp){l};
        s.(oname) = subplot(2,4,4+l);
        imagesc(images.plot.Hxbnds,images.plot.Hybnds,plotimg);
        title(oname,'FontSize',16); axis('equal','off'); hold on
        if axsl >= stHa && axsl <= fnHa %check if current slice is within 13C pH slab
            % Plot voxels for selected ROI and overlay
            rname = roinames{actroi};
            v = imagesc(images.plot.Pxbnds*pixratP,images.plot.Pybnds*pixratP,...
                images.(oname),lims.(actovlygrp)(l,:));
            if contains(oname,'pHe') %flip colorbar for pH
                colormap(s.(oname),flipud(jet)); colorbar;
            elseif ~strcmp(oname,'rescuedvoxels') %no colorbar for 
                %rescued voxel map
                colormap(s.(oname),'jet'); colorbar;
            end
            set(v,'AlphaData',rois.(rname).maskP*otp.*...
                (abs(images.(oname)) > 0)); 
            if ovlyflg && sum(sum(dispvox)) > 0 %outline selected voxels
                [olvoxy,olvoxx] = find(dispvox);
                xl = images.plot.Pxbnds(1) + images.CSI.ps(1) * (min(olvoxx)-1);
                xu = images.plot.Pxbnds(1) + images.CSI.ps(1) * max(olvoxx);
                yl = images.plot.Pybnds(1) + images.CSI.ps(1) * (min(olvoxy) - 1);
                yu = images.plot.Pybnds(1) + images.CSI.ps(1) * max(olvoxy);
                plot([xl xu xu xl xl],[yl yl yu yu yl],'g-','LineWidth',2);
            end
        end
        % Plot and label scale bar        
        plotaxes = [images.plot.Hxbnds,images.plot.Hybnds] .* ones(1,4) ...
            * (1-zoomfac);
        barx = [plotaxes(2)*0.9-4 plotaxes(2)*0.9];
        bary = [plotaxes(4) plotaxes(4)] * 0.9;
        plot(barx,bary,'w-','LineWidth',4)
        text(mean(barx),bary(1)*.8/.9,'4 mm','FontSize',12,'Color','w',...
            'HorizontalAlignment','center');
        axis(plotaxes); %adjust zoom
        hold off
    end
end

% plotCorOvly: plots coronal 1H image (+ overlay voxels, if within 13C pH 
% slab)
%
function plotCorOvly()
    % First, convert image to grayscale (necessary for overlays)
    plotimg = squeeze(images.cor.img(:,:,corsl));
    hlims = [min(min(plotimg)) max(max(plotimg))];
    plotimg = repmat(mat2gray(plotimg,hlims),[1 1 3]);
    for l = 1:2 %cycle through all overlays specified in active ovlynames group
        oname = ovlynames.(actovlygrp){l};
        s.(oname) = subplot(2,4,4+l); 
        imagesc(images.plot.Kxbnds,images.plot.Kzbnds,plotimg); 
        title(oname,'FontSize',16); axis('equal','off'); hold on
        if corsl >= stHc && corsl <= fnHc %check if current slice is within 13C pH FOV
            % Plot voxels for selected ROI and overlay
            rname = roinames{actroi};
            v = imagesc(images.plot.Pxbnds*pixratP,images.plot.Pzbnds*pixratz,...
                images.(oname)(yindP,:),lims.(actovlygrp)(l,:));
            if contains(oname,'pHe') %flip colorbar for pH
                colormap(s.(oname),flipud(jet)); colorbar;
            else
                colormap(s.(oname),'jet'); colorbar;
            end
            set(v,'AlphaData',rois.(rname).maskP(yindP,:)*otp.*...
                (abs(images.(oname)(yindP,:)) > 0));
        end
        % Plot and label scale bar        
        plotaxes = [images.plot.Kxbnds,images.plot.Kzbnds] .* ones(1,4) ...
            * (1-zoomfac);
        barx = [plotaxes(2)*0.9-4 plotaxes(2)*0.9];
        bary = [plotaxes(4) plotaxes(4)] * 0.9;
        plot(barx,bary,'w-','LineWidth',4)
        text(mean(barx),bary(1)*.8/.9,'4 mm','FontSize',12,'Color','w',...
            'HorizontalAlignment','center');
        axis(plotaxes); %adjust zoom
        hold off
    end
end


%% INTERNAL FUNCTIONS - GUI INTERACTION
%
% axSelect: Plot selected axial slice
%
function axSelect(source,~)
    axsl = round(source.Value);
    plotAxImg;
    if ovlyflg
        plotAxOvly;
    end
end

% corSelect: Plot selected coronal slice
%
function corSelect(source,~)
    corsl = round(source.Value);
    yindP = round(images.CSI.npx(2) / 2 + 0.5 - (images.cor.off(2) ...
        - images.CSI.off(2) + (corsl - images.cor.np(3) / 2 - 0.5) ...
        * images.cor.ps(3)) / images.CSI.ps(2));  
        %row index of pH data corresponding with active coronal slice
    plotCorImg;
    if ovlyflg
        plotCorOvly;
    end
end

% setAx: Sets active axes to axial figure
%
function setAx(~,~)
    curax = a; 
    plotAxImg; plotCorImg; %this will break out of the active imfreehand instance
end

% setCor: Sets active axes to coronal figure
%
function setCor(~,~) 
    curax = c;
    plotAxImg; plotCorImg; %this will break out of the active imfreehand instance
end

% setBounds: Toggle whether new integral bounds are defined for denoised
% data once data processing figure (specfig) is done with
%
function setBounds(source,~)
if source.Value == 1
    boundflg = true;
else
    boundflg = false;
end
end

% setAlign: Toggle whether to automatically align peaks 
%
function setAlign(source,~)
if source.Value == 1
    alignflg = true;
    %then need to proc CSI w/o shifts, calculate shiftmap, then apply to 
    %data
    procCSI({num2str(oldnp(1)),num2str(oldnp(2)),num2str(images.CSI.apod), ...
        num2str(diffvox(1)),num2str(diffvox(2))},true,false);
    calcMaps(false,true,true);
else
    alignflg = false;
end
procCSI({num2str(oldnp(1)),num2str(oldnp(2)),num2str(images.CSI.apod), ...
    num2str(diffvox(1)),num2str(diffvox(2))},true,alignflg);    
calcMaps;
end

% nameROI: sets name for new ROI once created
%
function nameROI(source,~)
    newname = source.String;
end

% newROI: initializes new ROI, 
%
function newROI(~,~)
    nROI = nROI + 1;
    nROIctr = nROIctr + 1;
    actroi = nROI;
    roinames{actroi} = newname;
    name = roinames{actroi};
    set(rl,'String',roinames); %updates list of ROI names for active ROI selection
    newname = ['ROI' num2str(nROIctr+1)];
    set(nr,'String',newname); %updates text in ROI naming box
    rois.(name).maskP = false(images.CSI.npx(1));
    editROI;
end

% editROI: llows user to select voxels for ROI until user clicks button to 
% accept ROI edits
%
function editROI(~,~)
    name = roinames{actroi};
    endflg = false;
    while 1 %this loop is only broken once the "Accept ROI edits" button is clicked
        if endflg
            break
        end
        h = drawfreehand(curax,'Closed',false); %draw freehand to select voxels
        % ID and highlight selected pH voxels on figures
        if isvalid(h)
            coords = h.Position;
            delete(h);
            mask = false(images.CSI.npx(1));
            for l = 1:size(coords,1) %create mask to ID pH voxels that freehand went through
                for m = 1:images.CSI.npx(1)
                    xl = images.plot.Pxbnds(1) + images.CSI.ps(1) * (m - 1);
                    xu = images.plot.Pxbnds(1) + images.CSI.ps(1) * m;
                    for n = 1:images.CSI.npx(2)
                        yl = images.plot.Pybnds(1) + images.CSI.ps(1) * (n - 1);
                        yu = images.plot.Pybnds(1) + images.CSI.ps(1) * n;
                        if isequal(curax,a) %mask is created differently if drawn on axial vs coronal
                            if coords(l,1) > xl && coords(l,1) < xu && coords(l,2) > yl ...
                                    && coords(l,2) < yu
                                mask(n,m) = 1;
                            end
                        else
                            if coords(l,1) > xl && coords(l,1) < xu && coords(l,2) > images.plot.Pzbnds(1) ...
                                    && coords(l,2) < images.plot.Pzbnds(2)
                                mask(yindP,m) = 1;
                            end
                        end
                    end
                end
            end
            rois.(name).maskP = logical(abs(rois.(name).maskP - mask)); 
            %reverse all mask values for selected voxels
            plotAxImg; plotCorImg;
        end
    end
end

% endROIEdits: Finishes any edits currently being performed on ROIs
%
function endROIEdits(~,~)
    endflg = true;
    plotAxImg; plotCorImg; %this will break out of the active imfreehand instance
end

% delROI: Delete active ROI, replot
%
function delROI(~,~)
    endROIEdits;
    nROI = nROI - 1;
    delname = roinames{actroi};
    rois = rmfield(rois,delname);
    roinames{actroi} = [];
    roinames = roinames(~cellfun('isempty',roinames));
    actroi = actroi - 1;
    set(rl,'String',roinames); %updates list of ROI names for active ROI selection
    plotAxImg; plotCorImg;
end

% selIntegralBounds: Toggle integral bound display on spectra on/off
%
function selIntegralBounds(~,~)
    bndplotflg = ~bndplotflg;
    if bndplotflg %update button string based on value
        set(ib,'String','Turn integral bounds off');
    else
        set(ib,'String','Turn integral bounds on');
    end
    plotSpecVox(bndplotflg);
end

% selGrid: Toggle voxel grid on/off
%
function selGrid(~,~)
    gridflg = ~gridflg;
    if gridflg %update button string based on value
        set(gb,'String','Turn grid off');
    else
        set(gb,'String','Turn grid on');
    end
    plotAxImg; plotCorImg;
    if ovlyflg
        plotAxOvly;
    end
end

% selVoxOvly: Toggle selected voxel overlay on/off
%
function selVoxOvly(~,~)
    voxflg = ~voxflg;
    if voxflg %update button string based on value
        set(vob,'String','Turn voxel overlay off');
    else
        set(vob,'String','Turn voxel overlay on');
    end
    plotAxImg; plotCorImg;
    if ovlyflg
        plotAxOvly;
    end
end

% setVoxTransparency: Adjust transparency of selected voxel overlay
%
function setVoxTransparency(source,~)
    vstp = source.Value;
    plotAxImg; plotCorImg;
    if ovlyflg
        plotAxOvly;
    end
end

% setOvlyTransparency: Adjust transparency of all 13C overlays
%
function setOvlyTransparency(source,~)
    otp = source.Value;
    plotAxImg; 
    if ~specflg
        plotCorImg;
    end
    if ovlyflg
        plotAxOvly;
    end
end

% setZoom: Zoom in on all images in final figure
%
function setZoom(source,~)
    zoomfac = source.Value;
    plotAxImg; plotCorImg;
    if ovlyflg
        plotAxOvly;
    end
end

% resetZoom: Reset zoom to original amount
%
function resetZoom(~,~)
    zoomfac = 0;
    plotAxImg; plotCorImg;
    if ovlyflg
        plotAxOvly;
    end
end

% selROI: Set active ROI for editing/display purposes
%
function selROI(source,~)
    actroi = source.Value;
    plotAxImg; plotCorImg;
    if ovlyflg
        plotAxOvly;
    end
end

% selOvly: Set active overlay for display purposes (for spectral figure
% only, prior to pH map calculation)
%
function selOvly(source,~)
ovlyname = ovlynamesproc{source.Value};
plotAxImg;
end

% selOvlyGrp: Set active overlay group (1 = urea+LacPyrRatio+pHe, 
% 2 = HP001+Pyr+Lac) for display purposes
%
function selOvlyGrp(source,~)
actovlygrp = ovlygrps{source.Value};
plotAxOvly;
end

% finishFig: Finish all figure stuff, continue function
%
function finishFig(~,~)
finishflg = true;
specflg = false; %end specfig stuff
if exist('specfig','var')
    if isvalid(specfig)
        close(specfig)
    end
end
end

% finishROIs: Finish all ROI editing, continue function
%
function finishROIs(~,~)
    endROIEdits;
    close all
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
dataSize = size(image);
% options for phase correction search (if used)
% bruteForce = 0;
simplexZeroOrder = 1;
% simplexZeroAndFirstOrder = 2;

params.phaseSensitive = false;
params.phaseCorrectTimePoints = false; 
params.nonNegativePenalty = true;
params.searchMethod = simplexZeroOrder;
params.baselineCorrect = false;
params.specPoints = dataSize(1);
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
    y = sv{ii};
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
%% OLD CODE
% When user initially was able to specify a single set of integral bounds
% to use for both original and denoised datasets:
% 
%     % Prompt user to specify integral bounds for 2D CSI data (or leave
%     % blank if user wants to draw manually on spectra)
%     %
%     prompt = {[spnames{1} ': integral start (ppm):'],...
%         [spnames{1} ': integral end (ppm):'],...
%         [spnames{2} ': integral start (ppm):'],...
%         [spnames{2} ': integral end (ppm):']};
%     dlg_title = 'Integration (leave any value blank to define manually)';
%     num_lines = 1;
%     def = {'' , '' , '' , ''};    
%     answer3 = inputdlg(prompt , dlg_title , num_lines , def);
%
%     if sum(strcmp(answer3,'')) == 0 %detect if all integral bounds were defined
%         if str2double(answer3{1}) < str2double(answer3{2}) %check which is 
%             %greater, which is lesser
%             images.spbin.(spnames{1}) = find(images.ppm > str2double(answer3{1}) & images.ppm < ...
%                 str2double(answer3{2}));
%         else
%             images.spbin.(spnames{1}) = find(images.ppm > str2double(answer3{2}) & images.ppm < ...
%                 str2double(answer3{1}));
%         end
%         if str2double(answer3{3}) < str2double(answer3{4})
%             images.spbin.(spnames{2}) = find(images.ppm > str2double(answer3{3}) & images.ppm < ...
%                 str2double(answer3{4}));
%         else
%             images.spbin.(spnames{2}) = find(images.ppm > str2double(answer3{4}) & images.ppm < ...
%                 str2double(answer3{3}));
%         end
%     end