%
% Batch parcellation resampling and source-reconstruction 
% for MEG HCP data analyis
%
% Alessandra Griffa
% University of Geneva, Switzerland
% April 2020
%

clearvars
close all
clc


% SET PATHS ---------------------------------------------------------------
% Update PATH - Include wb_command
% Download workbench: https://www.humanconnectome.org/software/get-connectome-workbench
setenv('PATH', [getenv('PATH') ':/Users/alli/Documents/CODE/hcp/workbench/bin_macosx64']);

% Addpath - megconnectom HCP functions, fieldtrip and cifti-matlab
% Download megconnectome, fieldtrip: https://www.humanconnectome.org/software/hcp-meg-pipelines
% Download cifti-matlab: https://github.com/Washington-University/cifti-matlab 
addpath(genpath('/Users/alli/Documents/CODE/hcp/megconnectome-3.0'));
addpath('/Users/alli/Documents/CODE/hcp/fieldtrip-r10442');
ft_defaults
addpath(genpath('/Users/alli/Documents/CODE/hcp/cifti-matlab-master'));

% Addpath script folder
addpath('./');

% Paths 4k-vertex Sphere surfaces, distributed with megconnectome
sphereR = '/Users/alli/Documents/CODE/hcp/megconnectome-3.0/template/Sphere.4k.R.surf.gii';
sphereL = '/Users/alli/Documents/CODE/hcp/megconnectome-3.0/template/Sphere.4k.L.surf.gii';


% SET MEG DATA PATH and OUTPUT DIR ----------------------------------------
data_path = '/Users/alli/Desktop/example_data/DATA_HCP_MEG';    % input dir
mat_path = '/Users/alli/Desktop/example_data/DATA_MAT';         % output dir
subj_list = dir(data_path);
% remove all files (isdir property is 0)
subj_list = subj_list([subj_list(:).isdir]==1);
% remove '.' and '..'
subj_list = subj_list(~ismember({subj_list(:).name},{'.','..'}));
% remove 'MAT' dir
subj_list = subj_list(~ismember({subj_list(:).name},{'MAT'}));

% Number of subjects
ns = length(subj_list);
disp([num2str(ns) ' subjects found in input data path ']);


% PARAMETER SETTING -------------------------------------------------------
% Do resample DK80 parcellation from 164k to 8k surface resolution?
DO_RESAMPLE = 0;
% Do perform source reconstruction?
DO_SOURCERECONSTRUCTION = 1;
% Set which Restin-preproc MEG sessions you want to reconstruct
session = [3, 4, 5];



%% Resample DK80 parcellation from 164k to 8k surface resolution
if DO_RESAMPLE
    
    % Loop over subjects
    for s = 1:ns

        subj = subj_list(s);

        disp([' > SURFACE RESMAPLING subject ' num2str(s) ' of ' num2str(ns) ' (id ' subj.name ')']);

        % -------- RESAMPLED PARCELLATION SAVED IN XXX_MEG_anatomy --------
        out_path = fullfile(data_path, subj.name, strcat(subj.name,'_MEG_anatomy'));

        % CHECK IF OUTPUT EXISTS
        if isfile(fullfile(out_path, strcat(subj.name, '.aparc.a2009s.8k_fs_LR.dlabel.nii')))
            disp('   output exists - subject skipped');
            continue
        end

        % ----------------- SET HERE PARCELLATION -------------------------
        this_parc_dlabel = fullfile(data_path, subj.name, strcat(subj.name,'_3T_Structural_preproc'), 'MNINonLinear', ...
                            strcat(subj.name,'.aparc.a2009s.164k_fs_LR.dlabel.nii'));

        this_L_sphere_surf = fullfile(data_path, subj.name, strcat(subj.name,'_3T_Structural_preproc'), 'MNINonLinear', ...
                            strcat(subj.name,'.L.sphere.164k_fs_LR.surf.gii')); 
        this_R_sphere_surf = fullfile(data_path, subj.name, strcat(subj.name,'_3T_Structural_preproc'), 'MNINonLinear', ...
                            strcat(subj.name,'.R.sphere.164k_fs_LR.surf.gii'));

        cmd = ['wb_command -cifti-separate ' this_parc_dlabel ' COLUMN  -label CORTEX_LEFT ' fullfile(out_path,'TMPOLD_LEFT.aparc.164k.label.gii')];
        status = system(cmd);
        cmd = ['wb_command -cifti-separate ' this_parc_dlabel ' COLUMN  -label CORTEX_RIGHT ' fullfile(out_path,'TMPOLD_RIGHT.aparc.164k.label.gii')];
        status = system(cmd);

        cmd = ['wb_command -label-resample ' fullfile(out_path,'TMPOLD_LEFT.aparc.164k.label.gii') ' ' ...
            this_L_sphere_surf ' ' sphereL ' BARYCENTRIC ' fullfile(out_path,'TMPNEW_LEFT.aparc.164k.label.gii')];
        status = system(cmd);
        cmd = ['wb_command -label-resample ' fullfile(out_path,'TMPOLD_RIGHT.aparc.164k.label.gii') ' ' ...
            this_R_sphere_surf ' ' sphereR ' BARYCENTRIC ' fullfile(out_path,'TMPNEW_RIGHT.aparc.164k.label.gii')];
        status = system(cmd);

        % -------------- SET HERE OUTPUT PARCELLATION FILE --------------------
        cmd = ['wb_command -cifti-create-label ' ...
            fullfile(out_path, strcat(subj.name, '.aparc.a2009s.8k_fs_LR.dlabel.nii')) ...
            ' -left-label ' fullfile(out_path,'TMPNEW_LEFT.aparc.164k.label.gii') ...
            ' -right-label ' fullfile(out_path,'TMPNEW_RIGHT.aparc.164k.label.gii')];
        status = system(cmd);

        % Clean temp files
        delete(fullfile(out_path,'TMPNEW_LEFT.aparc.164k.label.gii'));
        delete(fullfile(out_path,'TMPNEW_RIGHT.aparc.164k.label.gii'));
        delete(fullfile(out_path,'TMPOLD_LEFT.aparc.164k.label.gii'));
        delete(fullfile(out_path,'TMPOLD_RIGHT.aparc.164k.label.gii'));

    end
end



%% Source reconstruction
if DO_SOURCERECONSTRUCTION

    % Loop over subejcts
    for s =1:ns

        subj = subj_list(s);

        disp([' > SOURCE RECONSTRUCTION subject ' num2str(s) ' of ' num2str(ns) ' (id ' subj.name ')']);

        % Check if headmodel and sourcemodel files exist (they are provided by HCP)
        headmodelfile = fullfile(data_path, subj.name, strcat(subj.name,'_MEG_anatomy'), strcat(subj.name,'_MEG_anatomy_headmodel.mat'));
        sourcemodelfile = fullfile(data_path, subj.name, strcat(subj.name,'_MEG_anatomy'), strcat(subj.name,'_MEG_anatomy_sourcemodel_2d.mat'));
        if ~isfile(headmodelfile) || ~isfile(sourcemodelfile)
            disp(['  NO MODEL DATA for subject ' subj.name ]);
            continue
        end

        % Check if MEG data exists for session (preprocessed sensor-level MEG data provided by HCP)
        new_session = [];
        for ses = 1:length(session)
            megfile = fullfile(data_path, subj.name, strcat(subj.name,'_MEG_Restin_preproc'), 'rmegpreproc', strcat(subj.name,'_MEG_',num2str(session(ses)),'-Restin_rmegpreproc.mat'));
            outfile = dir(fullfile(mat_path, strcat(subj.name,'_SourceRecon_',num2str(session(ses)),'-Restin_rmegpreproc_aparc-a2009s_*')));
            if ~isfile(megfile)
                disp(['  NO DATA for subject ' subj.name ', Restin-preproc session ' num2str(session(ses))]);
            elseif ~isempty(outfile)
                disp(['  SOURCE RECONSTRUCTION EXISTIS for subject ' subj.name ', Restin-preproc session ' num2str(session(ses))]);
            else
                new_session = [new_session, session(ses)];
            end
        end
            
        if ~isempty(new_session) %&& isempty(outfile)
            % Call function
            fcn_source_recon_from_Restin_CorrectionDipole_Centroid(data_path, subj.name, new_session, mat_path);
        else
            disp(['  > > > subject ' subj.name ' skipped ']);
        end

    end
    
end


