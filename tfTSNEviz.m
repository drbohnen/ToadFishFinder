imds = imageDatastore('training', 'LabelSource', 'foldernames', 'IncludeSubfolders',true);

net = resnet50();
% get training features 
featureLayer = 'fc1000';
trainingFeatures = activations(net,imds, featureLayer, ...
    'MiniBatchSize', 32, 'OutputAs', 'columns');


% Get training labels from the trainingSet
trainingLabels = imds.Labels;

Y = tsne(trainingFeatures'); figure; gscatter(Y(:,1),Y(:,2),trainingLabels)