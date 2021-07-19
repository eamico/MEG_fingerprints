% meg_id_icc() code is associated with the manuscript titled 
% "Exploring MEG brain fingerprints: Evaluation, pitfalls, and
% interpretations" by Sareen et al., published in NeuroImage 
% (Volume 240, 15 October 2021, 118331). This example code can be 
% used to reproduce certain results from the manuscript. 
%
% Code Authors: Ekansh Sareen, Alessandra Griffa, Enrico Amico
% version 1.0 (15 July, 2021)
%
% Please cite us: 
% Ekansh Sareen, Sélima Zahar, Dimitri Van De Ville, Anubha Gupta, Alessandra Griffa, Enrico Amico,
% Exploring MEG brain fingerprints: Evaluation, pitfalls, and interpretations, NeuroImage,
% Volume 240, 2021, 118331, ISSN 1053-8119,https://doi.org/10.1016/j.neuroimage.2021.118331.
%
%% 0. Add relevant functions in the path

addpath(genpath('functions'));

%% [Part 1.1] Computing Identifiability Matrix
% Computing the subjectxsubject Identifiability matrix

% Load the session 1 and session 2 connectomes from the 'FCMethods' folder for
% 3 FC measures (AECc, PLM, wPLI). 

load(fullfile('FCMethods','AECc','tensor_sub_connmat_sesh1.mat'));
load(fullfile('FCMethods','AECc','tensor_sub_connmat_sesh2.mat'));
load(fullfile('FCMethods','PLM','tensor_sub_connmat_sesh1.mat'));
load(fullfile('FCMethods','PLM','tensor_sub_connmat_sesh2.mat'));
load(fullfile('FCMethods','wPLI','tensor_sub_connmat_sesh1.mat'));
load(fullfile('FCMethods','wPLI','tensor_sub_connmat_sesh2.mat'));

% Connectomes: 148x148x5x20 (148 brain regions, 5 frequency bands, 20 subjects)

bands=5;    % No. of frequency bands
sub=20;     % No. of Subjects
% For AECc
tensor1=tensor_sub_connmat_sesh1_AECc;
tensor2=tensor_sub_connmat_sesh2_AECc;
[~,~,ID_mat_AECc]= ID_construct(tensor1,tensor2,bands,sub);

% For PLM
tensor3=tensor_sub_connmat_sesh1_PLM;
tensor4=tensor_sub_connmat_sesh2_PLM;
[~,~,ID_mat_PLM]= ID_construct(tensor3,tensor4,bands,sub);

% For wPLI
tensor5=tensor_sub_connmat_sesh1_wPLI;
tensor6=tensor_sub_connmat_sesh2_wPLI;
[~,~,ID_mat_wPLI]= ID_construct(tensor5,tensor6,bands,sub);

% Visualization 
plot_graph(ID_mat_AECc,bands,sub,'AECc')
plot_graph(ID_mat_PLM,bands,sub,'PLM')
plot_graph(ID_mat_wPLI,bands,sub,'wPLI')

%% [Part 1.2] Histograms: Idiff and SR. Please note that SR values might vary depending on the sample size.

list_aecc=id_params_mod(ID_mat_AECc,bands,sub);
list_plm=id_params_mod(ID_mat_PLM,bands,sub);
list_wpli=id_params_mod(ID_mat_wPLI,bands,sub);

Idiff_aecc=list_aecc(3,:); Idiff_plm=list_plm(3,:); Idiff_wpli=list_wpli(3,:);
sr_aecc=list_aecc(4,:); sr_plm=list_plm(4,:); sr_wpli=list_wpli(4,:);

% Histogram of Idiff values for 3 FC measures and 5 frequency bands
figure;
h=[Idiff_aecc;Idiff_plm;Idiff_wpli];
g=bar(h);
g(1).FaceColor=[0.313 1 0];
g(2).FaceColor=[0.875 0 1];
g(3).FaceColor=[1 0.938 0];
g(4).FaceColor=[1 0 0];
g(5).FaceColor=[0 0.813 1];
title('Idiff Scores')
name={'AECc';'PLM';'wPLI'};
set(gca,'xticklabel',name);
set(g, {'DisplayName'}, {'Delta','Theta','Alpha', 'Beta', 'Gamma'}')
legend() 

% Histogram of SR values for 3 FC measures and 5 frequency bands
figure;
f=[sr_aecc;sr_plm;sr_wpli];
j=bar(f);
j(1).FaceColor=[0.313 1 0];
j(2).FaceColor=[0.875 0 1];
j(3).FaceColor=[1 0.938 0];
j(4).FaceColor=[1 0 0];
j(5).FaceColor=[0 0.813 1];
title('SR Scores (in %)')
name={'AECc';'PLM';'wPLI'};
set(gca,'xticklabel',name);
set(j, {'DisplayName'}, {'Delta','Theta','Alpha', 'Beta', 'Gamma'}')
legend() 

%% [Part 2] ICC Visualization
% ICC matrices for 3 FC measures and 2 frequency bands
icc_gen(tensor1,tensor2,'AECc');
icc_gen(tensor3,tensor4,'PLM');
icc_gen(tensor5,tensor6,'wPLI');
