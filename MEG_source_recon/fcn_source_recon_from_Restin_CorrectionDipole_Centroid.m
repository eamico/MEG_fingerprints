function [] = fcn_source_recon_from_Restin_CorrectionDipole_Centroid(cdbDir, subjectid, session, outDir)
%
% This function performs the source reconstruction of sensor-level MEG data
% from the HCP database, and following the script by Georgios Michalareas 
% (https://www.mail-archive.com/hcp-users@humanconnectome.org/msg06067.html)
% and the FieldTrip tutorial
% http://www.fieldtriptoolbox.org/tutorial/networkanalysis/.
% The source-resonstruciton is performed using the realistic single shell
% volume conduciton model (Nolte, 2003) and LCMV beamformer, with respect 
% to 8K vertices on the cortical surface. 
% Next, the time series of the vertices belonging to individual fs_a2009s 
% cortical parcels are averaged to provide a source-reconstructed time 
% series for each brain region.
% The region-level (source) time series are save in an output mat file.
%
% INPUTS
% - cdbDir: input data path (input data were downloaded from the HCP 
%           database) [string]
% - subjectid: subject id (e.g., '100307') [string]
% - session: rs-MEG reconrding section. Typically, each subject has 3
%            recording sessions, numbered 3 / 4 / 5  [int or inr array]
% - outDir: output directory where the outputs of this function will be
%           saved [string]
%
% Alessandra Griffa, alessandra.griffa@gmail.com
% EPFL, UNIGE
% March 2020
%


% Setting INSPECT=0 avoids plotting
INSPECT = 0;



%% LOAD DATA
% - sourcemodel2d
% - headmodel
% - data            channel-level time series
% -------------------------------------------------------------------------


% Load SOURCEMODEL2D (source coordinates ~ 8k vertices of cortical mesh)
% -------------------------------------------------------------------------
% 'bti' coordinate system: see http://www.fieldtriptoolbox.org/faq/how_are_the_different_head_and_mri_coordinate_systems_defined/
anatomydir = fullfile(cdbDir, subjectid, strcat(subjectid,'_MEG_anatomy'));
fileGrid2D = fullfile(anatomydir, strcat(subjectid,'_MEG_anatomy_sourcemodel_2d.mat'));
hcp_read_matlab(fileGrid2D,'sourcemodel2d');    % HCP_READ_MATLAB reads one or multiple MATLAB variable from a *.mat file.
gridd = ft_convert_units(sourcemodel2d,'cm');   % FT_CONVERT_UNITS changes the geometrical dimension to the specified SI unit


% SOURCEMODEL2DROI: generate a new source model considering the centroids 
%                   of the anatomical regions 
% -------------------------------------------------------------------------
% Load parcellation 
this_parc = fullfile(cdbDir, subjectid, strcat(subjectid,'_MEG_anatomy'), ...
            strcat(subjectid,'.aparc.a2009s.8k_fs_LR.dlabel.nii'));
parc = ciftiopen(this_parc);
parc_id = sort(unique(parc.cdata));
parc_id = parc_id(parc_id ~= 0);
% Number of regions in current parcellation
nrois = length(parc_id);
% Build new sourcemodel sstructure
sourcemodel2dROI.coordsys = sourcemodel2d.coordsys;
sourcemodel2dROI.unit = sourcemodel2d.unit;
sourcemodel2dROI.brainstructurelabel = sourcemodel2d.brainstructurelabel;
sourcemodel2dROI.pos = zeros(nrois,3);
sourcemodel2dROI.brainstructure = zeros(nrois,1);
T = parc.diminfo{2}.maps.table;
Tlabels = extractfield(T,'name');
Tkeys = extractfield(T,'key');
for i = 1:nrois
    ii = find(parc.cdata == parc_id(i));
    temp = sourcemodel2d.pos(ii,:);
    sourcemodel2dROI.pos(i,:) = mean(temp);
    ii = find(Tkeys == i);
    temp = Tlabels(ii);
    if contains(temp,'L_')
        sourcemodel2dROI.brainstructure(i) = 1;
    else
        sourcemodel2dROI.brainstructure(i) = 2;
    end
end
griddROI = ft_convert_units(sourcemodel2dROI,'cm');   % FT_CONVERT_UNITS changes the geometrical dimension to the specified SI unit


if INSPECT
    figure
    plot3(gridd.pos(gridd.brainstructure==1, 1), ...
            gridd.pos(gridd.brainstructure==1, 2), ...
            gridd.pos(gridd.brainstructure==1, 3), 'ro'), hold on;
    plot3(gridd.pos(gridd.brainstructure==2, 1), ...
            gridd.pos(gridd.brainstructure==2, 2), ...
            gridd.pos(gridd.brainstructure==2, 3), 'co'), hold on;
    axis equal tight, grid on, set(gcf,'color','w'), set(gca,'fontsize',18), title('SOURCE MODEL');
    
    figure
    plot3(griddROI.pos(griddROI.brainstructure==1, 1), ...
            griddROI.pos(griddROI.brainstructure==1, 2), ...
            griddROI.pos(griddROI.brainstructure==1, 3), 'ro'), hold on;
    plot3(griddROI.pos(griddROI.brainstructure==2, 1), ...
            griddROI.pos(griddROI.brainstructure==2, 2), ...
            griddROI.pos(griddROI.brainstructure==2, 3), 'co'), hold on;
    axis equal tight, grid on, set(gcf,'color','w'), set(gca,'fontsize',18), title('SOURCE MODEL CENTROIDS');
end


% Load HEADMODEL
% 'bnd' contains the geometrical description of the head model
%       The bnd field describes a surface with vertices and triangles (in the bnd.pnt and bnd.tri fields)
%       as the geometrical description of the volume conductor.
% 'type' describes the method that was used to create the headmodel.
% 'unit' the unit of measurement of the geometrical data in the bnd field
% 'cfg' configuration of the function that was used to create the headmodel
% Note that this realistic single-shell model was created using brain surface from segmented mri.
% -------------------------------------------------------------------------
fileHeadVol = fullfile(anatomydir, strcat(subjectid,'_MEG_anatomy_headmodel.mat'));
hcp_read_matlab(fileHeadVol,'headmodel');
headmodel = ft_convert_units(headmodel,'cm');

if INSPECT
    figure, trimesh(headmodel.bnd.tri, headmodel.bnd.pnt(:,1), headmodel.bnd.pnt(:,2), headmodel.bnd.pnt(:,3));
    axis equal tight, grid on, set(gcf,'color','w'), set(gca,'fontsize',18), title('HEAD MODEL');
end


% Loop over sessions
for ses = 1:length(session)
    
    
    % LOAD CHANNEL-LEVEL TIME SERIES
    % -------------------------------------------------------------------------
    filename = fullfile(cdbDir, subjectid, strcat(subjectid,'_MEG_Restin_preproc'), 'rmegpreproc', strcat(subjectid,'_MEG_',num2str(session(ses)),'-Restin_rmegpreproc.mat'));
    hcp_read_matlab(filename,'data');
    grad = ft_convert_units(data.grad,'cm'); % assuming that subject has moved very little between the two scans

    if INSPECT
        figure, plot3(grad.coilpos(1:241,1),grad.coilpos(1:241,2),grad.coilpos(1:241,3),'b*');
        axis equal tight, grid on, set(gcf,'color','w'), set(gca,'fontsize',18), title('SENSOR POSITION');

        % Plot the source model, the head model and the sensor position in the same space
        figure
        plot3(gridd.pos(gridd.brainstructure==1, 1), ...
                gridd.pos(gridd.brainstructure==1, 2), ...
                gridd.pos(gridd.brainstructure==1, 3), 'ro'), hold on;
        plot3(gridd.pos(gridd.brainstructure==2, 1), ...
                gridd.pos(gridd.brainstructure==2, 2), ...
                gridd.pos(gridd.brainstructure==2, 3), 'co');
        trimesh(headmodel.bnd.tri, headmodel.bnd.pnt(:,1), headmodel.bnd.pnt(:,2), headmodel.bnd.pnt(:,3)), hold on;
        plot3(grad.coilpos(1:241,1),grad.coilpos(1:241,2),grad.coilpos(1:241,3),'b*');
        axis equal tight, grid on, set(gcf,'color','w'), set(gca,'fontsize',18), title('MEG SENSORS and HEAD MODEL');
    end




    %% CREATE LEAD FIELD MATRICES ~ gridAllLF
    % 8004 matrices > 241 x 3 (241: number of channels; 3: x-y-z directions of the dipole)
    % -------------------------------------------------------------------------
    % ( you only need to do this once per subject. Not everytime you perform
    % source localization ) 
    %
    % FT_CHANNELSELECTION makes a selection of MEG channel labels.
    % This function translates the user-specified list of channels into channel
    % labels as they occur in the data.

    allChansMEG = ft_channelselection({'MEG'},grad.label);


    % FT_PREPARE_LEADFIELD computes the forward model for many dipole locations
    % (gridd) and stores it for efficient inverse modelling
    cfg = [];
    cfg.grid = gridd;               % SOURCEMODEL ~ brain locations of the sources
    cfg.headmodel = headmodel;      % HEADMODEL ~ brain mesh
    cfg.grad = grad;                % SENSORS POSITINONS
    cfg.reducerank      = 2;        % (default = 3 for EEG, 2 for MEG)
    cfg.normalize       = 'yes' ;   % Normalise Leadfield: 'yes' for beamformer
    cfg.normalizeparam  = 1;        % depth normalization parameter (default = 0.5).
    cfg.feedback = 'no';
    cfg.channel = allChansMEG;
    cfg.inwardshift = -1;           % just to meake sure that the entire cortical sheet is in source space

    gridAllLF = ft_prepare_leadfield(cfg, data);    % gridAllLG.leadfield{1}: 241 x 3 matrix 
                                                    % ~ lead field spatial filter along the x-y-z dipole directions 
    gridAllLF.label = allChansMEG;  % CHANNELS LABELS
    
    
    % FT_PREPARE_LEADFIELD computes the forward model for BRAIN REGION
    % CENTROID DIPOLE LOCATIONS (griddROI) and stores it for efficient inverse 
    % modelling
    cfg = [];
    cfg.grid = griddROI;            % SOURCEMODEL ~ brain locations of the sources
    cfg.headmodel = headmodel;      % HEADMODEL ~ brain mesh
    cfg.grad = grad;                % SENSORS POSITINONS
    cfg.reducerank      = 2;        % (default = 3 for EEG, 2 for MEG)
    cfg.normalize       = 'yes' ;   % Normalise Leadfield: 'yes' for beamformer
    cfg.normalizeparam  = 1;        % depth normalization parameter (default = 0.5).
    cfg.feedback = 'no';
    cfg.channel = allChansMEG;
    cfg.inwardshift = -1;           % just to meake sure that the entire cortical sheet is in source space

    gridAllLFROI = ft_prepare_leadfield(cfg, data);    % gridAllLG.leadfield{1}: 241 x 3 matrix 
                                                    % ~ lead field spatial filter along the x-y-z dipole directions 
    gridAllLFROI.label = allChansMEG;  % CHANNELS LABELS


    if INSPECT

        CM = colormap(jet); close all
        val = rescale(gridAllLF.leadfield{1}(:,1)) * size(CM,1);
        figure;
        temp = ceil(val); temp(temp <= 0) = 1; temp(temp > size(CM,1)) = size(CM,1);
        scatter3(grad.coilpos(1:241,1),grad.coilpos(1:241,2),grad.coilpos(1:241,3), [], CM(temp,:), 'filled');
        hold on,
        scatter3(gridd.pos(1,1),gridd.pos(1,2),gridd.pos(1,3), 40, '*k');
        axis equal tight, grid on, set(gcf,'color','w'), set(gca,'fontsize',18), title('LEAD FIELD spatial filter for a single dipole (vertex reconstruction)');

        val = rescale(gridAllLFROI.leadfield{1}(:,1)) * size(CM,1);
        figure;
        temp = ceil(val); temp(temp <= 0) = 1; temp(temp > size(CM,1)) = size(CM,1);
        scatter3(grad.coilpos(1:241,1),grad.coilpos(1:241,2),grad.coilpos(1:241,3), [], CM(temp,:), 'filled');
        hold on,
        scatter3(griddROI.pos(1,1),griddROI.pos(1,2),griddROI.pos(1,3), 40, '*k');
        axis equal tight, grid on, set(gcf,'color','w'), set(gca,'fontsize',18), title('LEAD FIELD spatial filter for a single dipole (centroid reconstruction)');
        
    end




    %% COMPUTE SOURCE-LEVEL TIME SERIES (LCMV beamformer)
    % In addition to a forward model, the beamformer needs a sensor-level 
    % covariance matrix, or a cross-spectral density matrix (depending on the
    % choosen beamformer method; LCMV requires a data covariance matrix, DICS
    % requires a cross-spectral density matrix)

    % COMPUTE THE DATA COVARIANCE MATRIX
    % OPTION 1: average the time series across trials, then compute the
    %           covariance matrix over the averaged time series
    % OPTION 2: compute a covariance matrix for each trial and then average the
    %           covariance matrices (this is equivalent to computing the
    %           covariance matrix on the concatenated time series)
    OPTION = 2;

    if OPTION == 1
        % Compute average time series, so that the covariance matrix will be
        % computed on the 
        cfg = [];
        cfg.covariance = 'no';
        cfg.covariancewindow = 'all';
        datavgLH = ft_timelockanalysis(cfg, data);
    elseif OPTION == 2
        datavgLH = data;
    else
        disp('BAD OPTION FOR COVARIANCE MATRIX ESTIMATION');
    end

    % COMPUTE COVARIANCE MATRIX
    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = 'all';
    datavgLH = ft_timelockanalysis(cfg, datavgLH); 

    if INSPECT
        figure, imagesc(datavgLH.cov), colormap(jet), colorbar, xlabel('channels'), ylabel('channels');
        title('Channel-level data covariance matrix');
        set(gcf,'color','w'), set(gca,'fontsize',18), axis square;
    end


    %----- Set initial source localization settings
    % Linearly constrained minimum variance (LCMV) beamformers compute an estimate 
    % of source activity at each location through spatial filtering. The spatial 
    % data are linearly combined with weights (the spatial filter) chosen separately 
    % for each location to ensure that the strength of a dipolar source at that 
    % location is correctly estimated (assuming a perfect head model).
    srcCfg = [];                             
    srcCfg.headmodel = headmodel;      
    srcCfg.method = 'lcmv';
    srcCfg.lcmv.fixedori = 'yes';
    srcCfg.lcmv.feedback = 'text';
    srcCfg.lcmv.lambda = '20%';
    srcCfg.lcmv.keepfilter = 'yes';
    srcCfg.lcmv.keepmom = 'no';
    srcCfg.lcmv.projectnoise = 'no';
    srcCfg.keepleadfield = 'yes';
    srcCfgROI = srcCfg;

    % Select only channels of current data set in already computed leadfields 
    [indA,indB] = match_str(data.label,gridAllLF.label);
    inIndices = find(gridAllLF.inside);
    gridCase = gridAllLF;       % LEAD FIELD
    % Add grid (lead field) to source localization settings
    srcCfg.grid = gridCase;
    srcCfg.grid.label = gridCase.label(indB);
    
    % Select only channels of current data set in already computed
    % leadfields - ROI CENTROIDS
    [indA,indB] = match_str(data.label,gridAllLFROI.label);
    inIndices = find(gridAllLFROI.inside);
    gridCaseROI = gridAllLFROI;       % LEAD FIELD
    % Add grid (lead field) to source localization settings
    srcCfgROI.grid = gridCaseROI;
    srcCfgROI.grid.label = gridCase.label(indB);

    %----- Compute Inverse Solution from Covariance matrix ------------
    data.grad = grad;
    disp('beamforming starting...');
    sourceLocAll = ft_sourceanalysis(srcCfg, datavgLH);
    disp('...beamforming ended');
    
    disp('ROI beamforming starting...');
    sourceLocAllROI = ft_sourceanalysis(srcCfgROI, datavgLH);
    disp('...ROI beamforming ended');
    

    if INSPECT

        CM = colormap(jet); close all
        val = rescale(sourceLocAll.avg.pow) * size(CM,1);
        temp = ceil(val); temp(temp <= 0) = 1; temp(temp > size(CM,1)) = size(CM,1);
        figure;
        scatter3(sourceLocAll.pos(:,1),sourceLocAll.pos(:,2),sourceLocAll.pos(:,3), [], CM(temp,:), 'filled');
        axis equal tight, grid on, set(gcf,'color','w'), set(gca,'fontsize',18), title('Cortical power distribution');

        val = rescale(sourceLocAllROI.avg.pow) * size(CM,1);
        temp = ceil(val); temp(temp <= 0) = 1; temp(temp > size(CM,1)) = size(CM,1);
        figure;
        scatter3(sourceLocAllROI.pos(:,1),sourceLocAllROI.pos(:,2),sourceLocAllROI.pos(:,3), [], CM(temp,:), 'filled');
        axis equal tight, grid on, set(gcf,'color','w'), set(gca,'fontsize',18), title('Cortical power distribution - ROI centroids');

    end

    Nsensors = length(datavgLH.label);              disp(['>> Nsensors = ' num2str(Nsensors)]);
    Nsources = size(sourceLocAll.pos,1);            disp(['>> Nsources = ' num2str(Nsources)]);
    NsourcesROI = size(sourceLocAllROI.pos,1);      disp(['>> NsourcesROI = ' num2str(NsourcesROI)]);
    Ntimes = length(datavgLH.time);                 disp(['>> Ntimes = ' num2str(Ntimes)]);
    Ntrials = length(data.trial);                   disp(['>> Ntrials = ' num2str(Ntrials)]);


    % Initialize a matrix to keep all your single trial source leve data. 
    % This matrix can be large!!
    srcdata = nan(Nsources, Ntimes, Ntrials);
    srcdataROI = nan(NsourcesROI, Ntimes, Ntrials);

    % Put all spatial filters already computed in  a big matrix so that sensor
    % level data from each trial can be projected into all sources.
    allSpatialFilters = reshape([sourceLocAll.avg.filter{:}],Nsensors,Nsources)';
    allSpatialFiltersROI = reshape([sourceLocAllROI.avg.filter{:}],Nsensors,NsourcesROI)';
    if INSPECT
        figure, imagesc(allSpatialFilters), colormap(jet), colorbar, xlabel('channels'), ylabel('source locations');
        title([{'Spatial filters for back-projection'},{'(each row is a filter)'}]);
        set(gcf,'color','w'), set(gca,'fontsize',18);
    end


    %-------------------------------------------------------------------
    % Do the projection from sensor-level to source-level
    for iTrl = 1:Ntrials

        disp(['...back-projecting trial ' num2str(iTrl) ' of ' num2str(Ntrials)]);

        tmpdataTrial = data.trial{iTrl};
        srcdata(:,:,iTrl) = allSpatialFilters * tmpdataTrial;
        srcdataROI(:,:,iTrl) = allSpatialFiltersROI * tmpdataTrial;

    end
    %-------------------------------------------------------------------
    
    if INSPECT   
        test = srcdata(:,:,1);
        fc = corr(test');
        figure, imagesc(fc), colormap(jet), colorbar, axis square;
        title('FC Person correlation - trial 1');
        set(gcf,'color','w'); set(gca,'fontsize',18);
        
        test = srcdataROI(:,:,1);
        fc = corr(test');
        figure, imagesc(fc), colormap(jet), colorbar, axis square;
        title('FC Person correlation, ROI CENTROIDS - trial 1');
        set(gcf,'color','w'); set(gca,'fontsize',18);        
    end



    
    %% SAVE REGION-WISE TIME SERIES
        
    % Save data
    meg.subj = subjectid;
    meg.session = [num2str(session(ses)) '-Restin_rmegpreproc'];
    
    lastslash_pos = find(filename == '/', 1, 'last');
    meg.files.data = filename(lastslash_pos+1:end);
    lastslash_pos = find(fileHeadVol == '/', 1, 'last');
    meg.files.head = fileHeadVol(lastslash_pos+1:end);
    lastslash_pos = find(fileGrid2D == '/', 1, 'last');
    meg.files.source = fileGrid2D(lastslash_pos+1:end);
    meg.files.parc = parc.diminfo{2}.maps.name;
    
    meg.parc_id = parc_id;
    meg.parc_table = parc.diminfo{2}.maps.table;
    meg.fsample = data.fsample;
    meg.time = data.time;
    
    % Save ROI-CENTROID reconstruction
    % (instead of vertex-level data; other options have been tested)
    meg.ts_ROI = srcdataROI;

    time_stamp = datestr(now, 'yyyy-mm-dd');
    out_file = fullfile(outDir, strcat(meg.subj,'_SourceRecon_',num2str(session(ses)),'-Restin_rmegpreproc_aparc-a2009s_', time_stamp, '.mat'));
    save(out_file, 'meg'); 
   
end



