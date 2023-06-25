% This is the function use to train a new classifier 
% AUTHORS: 
% D. Bohnenstiehl (NCSU) 
% toadfish finder v.1 
% 29 Sept 22 

%% Load the scalogram images as a datastore, folder name  = labels 
imds = imageDatastore('TF_Training_v4', 'LabelSource', 'foldernames', 'IncludeSubfolders',true);
imds=shuffle(imds);  % mix the up 
tbl = countEachLabel(imds)

%% generate a sheet of examples 
bwhistle = find(imds.Labels == 'bwhistle');other = find(imds.Labels == 'other');
nP=table2array(tbl(1,2));
nO=table2array(tbl(2,2));
figure
for i=1:16 
subplot(4,4,i); imshow(readimage(imds,bwhistle(randi(nP,1)))); title('bwhistle example')
end
figure
for i=1:16 
subplot(4,4,i); imshow(readimage(imds,other(randi(nO,1)))); title('other example')
end

%% Use splitEachLabel method to trim the set.
minSetCount = min(tbl{:,2});
imds = splitEachLabel(imds, minSetCount, 'randomize');

%% Load pretrained network
net = resnet50();

%% Prepare Training and Test Image Sets
[trainingSet, testSet] = splitEachLabel(imds, 0.7, 'randomize');

%% Get Training Features 
featureLayer = 'fc1000';
trainingFeatures = activations(net,trainingSet, featureLayer, ...
    'MiniBatchSize', 32, 'OutputAs', 'columns');

%% Get training labels from the trainingSet
trainingLabels = trainingSet.Labels;

%% Define classifier options and trains the classifier.
classifier1 = fitcsvm(...
    trainingFeatures', ...
    trainingLabels, ...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 3, ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true, ...
    'Cost',[0 1;1 0],...
    'ClassNames', categorical({'other'; 'bwhistle'}))

%% estimate posterior probability from scores 
classifier = fitSVMPosterior(classifier1, 'KFold',10);  

%% False positive vs. true positive rate 
[~,scores1] = resubPredict(classifier);
[x,y,~,auc] = perfcurve(trainingLabels,scores1(:,2),'bwhistle');
auc
figure; 
plot(x,y)
xlabel('False positive rate') 
ylabel('True positive rate')
title('AUC')
ylim([0.975,1.0005])
xlim([-0.01,1]); grid on; 

%%  apply to training data 
[predictedLabels,TrainingScores] = predict(classifier, trainingFeatures, 'ObservationsIn', 'columns');

% Get the known labels
trainingLabels = trainingSet.Labels;

% Tabulate the results using a confusion matrix.
confMat = confusionmat(trainingLabels, predictedLabels);
disp('results for training data')

% Convert confusion matrix into percentage form
confMat = bsxfun(@rdivide,confMat,sum(confMat,2))


%% calculate k-folds cross model 
disp('K-fold Cross Validation Class Loss')
CVSVMModel = crossval(classifier1);
classLoss = kfoldLoss(CVSVMModel)

classifier2 = fitcsvm(...
    trainingFeatures', ...
    trainingLabels, ...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 3, ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true, ...
    'KFold',10,...
    'ClassNames', categorical({'other'; 'bwhistle'})); 

predictionkfold = kfoldPredict(classifier2);
confMat = confusionmat(trainingLabels, predictionkfold);
confMat = bsxfun(@rdivide,confMat,sum(confMat,2))

%%  apply to test data set 
% Extract test features using the CNN
testFeatures = activations(net, testSet, featureLayer, ...
    'MiniBatchSize', 32, 'OutputAs', 'columns');

% Pass CNN image features to trained classifier
[predictedLabels,predictedScores]= predict(classifier, testFeatures, 'ObservationsIn', 'columns');

% Get the known labels
testLabels = testSet.Labels;

% Tabulate the results using a confusion matrix.
confMat = confusionmat(testLabels, predictedLabels);

disp('results for test data')
% Convert confusion matrix into percentage form
confMat = bsxfun(@rdivide,confMat,sum(confMat,2))

figure; 
cm=confusionchart(testLabels, predictedLabels)
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
cm.Title = 'Test Data Confusion Matrix'

%% uncomment to save the classifier 
% save tfclassifer_vX.mat classifier 