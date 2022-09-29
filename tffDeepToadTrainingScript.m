% This is the function use to train a new classifier 
% AUTHORS: 
% D. Bohnenstiehl (NCSU) 
% toadfish finder v.1 
% 29 Sept 22 

%% Load the scalogram images as a datastore, folder name  = labels 
imds = imageDatastore('training', 'LabelSource', 'foldernames', 'IncludeSubfolders',true);
imds=shuffle(imds);  % mix the up 
tbl = countEachLabel(imds)

%% generate a sheet of examples 
boat = find(imds.Labels == 'boat');other = find(imds.Labels == 'other');

nP=table2array(tbl(1,2));
nO=table2array(tbl(2,2));
figure
for i=1:9 
subplot(3,3,i); imshow(readimage(imds,boat(randi(nP,1)))); title('boat example')
end

figure
for i=1:9 
subplot(3,3,i); imshow(readimage(imds,other(randi(nO,1)))); title('other example')
end


% Use splitEachLabel method to trim the set.
minSetCount = min(tbl{:,2});
imds = splitEachLabel(imds, minSetCount, 'randomize');

% Load pretrained network
net = resnet50();

%Prepare Training and Test Image Sets
[trainingSet, testSet] = splitEachLabel(imds, 0.6, 'randomize');


% Get Training Features 
featureLayer = 'fc1000';
trainingFeatures = activations(net,trainingSet, featureLayer, ...
    'MiniBatchSize', 32, 'OutputAs', 'columns');


% Get training labels from the trainingSet
trainingLabels = trainingSet.Labels;


% Define classifier options and trains the classifier.
classifier = fitcsvm(...
    trainingFeatures', ...
    trainingLabels, ...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 3, ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true, ...
    'ClassNames', categorical({'other'; 'boat'})); 


%%  apply to training data 
predictedLabels = predict(classifier, trainingFeatures, 'ObservationsIn', 'columns');

% Get the known labels
traiingLabels = trainingSet.Labels;

% Tabulate the results using a confusion matrix.
confMat = confusionmat(traiingLabels, predictedLabels);

disp('results for training data')
% Convert confusion matrix into percentage form
confMat = bsxfun(@rdivide,confMat,sum(confMat,2))


%%  apply to test data set 
% Extract test features using the CNN
testFeatures = activations(net, testSet, featureLayer, ...
    'MiniBatchSize', 32, 'OutputAs', 'columns');

% Pass CNN image features to trained classifier
predictedLabels = predict(classifier, testFeatures, 'ObservationsIn', 'columns');

% Get the known labels
testLabels = testSet.Labels;

% Tabulate the results using a confusion matrix.
confMat = confusionmat(testLabels, predictedLabels);

disp('results for test data')
% Convert confusion matrix into percentage form
confMat = bsxfun(@rdivide,confMat,sum(confMat,2))








