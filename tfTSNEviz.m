% This script reads a folder of labeled data with subfolder.
% Calculated ResNet-50 activations and displays a TSNE ordination  
% 
% ToadFishFinder v1.1 
% D. Bohnenstiehl 
% June 2023 
% 

imds = imageDatastore('TF_Training_v4', 'LabelSource', 'foldernames', 'IncludeSubfolders',true);
% Where 'TF_Training_v4' has subfolders 'other' and 'bwhistle' that hosted
% labeled spectrogram images. 

net = resnet50();
% get training features 
featureLayer = 'fc1000';
trainingFeatures = activations(net,imds, featureLayer, ...
    'MiniBatchSize', 32, 'OutputAs', 'columns');


% Get training labels from the trainingSet
trainingLabels = imds.Labels;

Y = tsne(trainingFeatures'); figure; gscatter(Y(:,1),Y(:,2),trainingLabels)