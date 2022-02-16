%
% Generate 4k parcellation dlabel 
%

clearvars
close all
clc



% Update PATH - Include wb_command
setenv('PATH', [getenv('PATH') ':/Users/alli/Documents/CODE/hcp/workbench/bin_macosx64']);



% PATH and SUBJECTS LIST
sphereR = '/Users/alli/Documents/CODE/hcp/megconnectome-3.0/template/Sphere.4k.R.surf.gii';
sphereL = '/Users/alli/Documents/CODE/hcp/megconnectome-3.0/template/Sphere.4k.L.surf.gii';

data_path = '/Users/alli/Desktop/2020_MEG_PCA/DATA';

subj_list = {'100307'};
ns = length(subj_list);



% Loop over subjects
for s = 1:ns
    
    subj = subj_list{s};
    
    disp([' > subject ' num2str(s) ' of ' num2str(s) ' (' subj ')']);
    
    % ----------------- SET HERE OUTPUT DIR -------------------------------
    out_path = fullfile(data_path, subj, strcat(subj,'_MEG_anatomy'), 'MEG', 'anatomy');
    
    % ----------------- SET HERE PARCELLATION -----------------------------
    this_parc_dlabel = fullfile(data_path, subj, strcat(subj,'_3T_Structural_preproc'), subj, 'MNINonLinear', ...
                        strcat(subj,'.aparc.a2009s.164k_fs_LR.dlabel.nii'));
    
    this_L_sphere_surf = fullfile(data_path, subj, strcat(subj,'_3T_Structural_preproc'), subj, 'MNINonLinear', ...
                        strcat(subj,'.L.sphere.164k_fs_LR.surf.gii')); 
    this_R_sphere_surf = fullfile(data_path, subj, strcat(subj,'_3T_Structural_preproc'), subj, 'MNINonLinear', ...
                        strcat(subj,'.R.sphere.164k_fs_LR.surf.gii'));
                    
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
        fullfile(out_path, strcat(subj, '.aparc.a2009s.8k_fs_LR.dlabel.nii')) ...
        ' -left-label ' fullfile(out_path,'TMPNEW_LEFT.aparc.164k.label.gii') ...
        ' -right-label ' fullfile(out_path,'TMPNEW_RIGHT.aparc.164k.label.gii')];
    status = system(cmd);
    
    % Clean temp files
    delete(fullfile(out_path,'TMPNEW_LEFT.aparc.164k.label.gii'));
    delete(fullfile(out_path,'TMPNEW_RIGHT.aparc.164k.label.gii'));
    delete(fullfile(out_path,'TMPOLD_LEFT.aparc.164k.label.gii'));
    delete(fullfile(out_path,'TMPOLD_RIGHT.aparc.164k.label.gii'));
    
end



