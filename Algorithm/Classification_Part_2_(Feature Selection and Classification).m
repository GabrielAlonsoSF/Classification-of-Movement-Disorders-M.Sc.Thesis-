clc
%{
features:
	standard deviation in time domain
	mean in time domain
	amplitude in time domain
	PSD peak
	dominant frequency
	second PSD peak 
	second dominant frequency
	cross-correlation

%}

% Since the gyroscope signals had a greater amplitude in all cases, only the gyroscope signals were used fot the classification

% creating a dataset for the features (observations X features)
all_feat = [gstd_X; gstd_Y; gstd_Z; gmean_X; gmean_Y; gmean_Z; gAmp_X; gAmp_Y; gAmp_Z; gAmp_SpX; gAmp_SpY; 
    gAmp_SpZ; gf_dom1X; gf_dom1Y; gf_dom1Z; gSp_2ndX; gSp_2ndY; gSp_2ndZ; gf_dom2X; gf_dom2Y; gf_dom2Z;
    gAmp_CRXY; gAmp_CRXZ; gAmp_CRYZ; g_tcorr_XY; g_tcorr_XZ; g_tcorr_YZ]';

% reading the correspondent class of the observations
c_012 = csvread('class_patients.csv');

% dataset joining the features and classes 
Dataset = [all_feat c_012];

% data balancing 
class_names = unique(c_012);
PD_cases = sum(c_012==0);
ES_cases = sum(c_012==1);
DS_cases = sum(c_012==2);

Dataset_v2 = sort(Dataset);
PD_Mat = Dataset_v2(1:PD_cases,:);
ES_Mat = Dataset_v2(PD_cases+1:PD_cases+ES_cases,:);
DS_Mat = Dataset_v2(PD_cases+ES_cases+1:end,:);
rng(10)
PD_perm = randperm(PD_cases);
rng(10)
ES_perm = randperm(ES_cases);
PD_Mat(PD_perm(1:PD_cases-DS_cases),:) = []; 
ES_Mat(ES_perm(1:ES_cases-DS_cases),:) = [];
all_feat = [PD_Mat; ES_Mat; DS_Mat];

% training and test division

s_feat = size(all_feat);
numObs = length(all_feat);
rng(10)
perm = randperm(numObs);
for m = 1:s_feat(2)
    for i=1:numObs
        all_feat_2(i,m) = all_feat(perm(i),m); %shuffling data for training and test division
    end
end

TrainPart = floor(numObs*0.7); % 70% of the observations
trn = all_feat_2(1:TrainPart,1:27); % training set
tst = all_feat_2(TrainPart+1:end,1:27); % test set 
trn_c = all_feat_2(1:TrainPart,28); % classes of training set
tst_c = all_feat_2(TrainPart+1:end,28); % classes of test set

% Feature selection using Principal Component Analysis (PCA)
 
all_sz = size(all_feat_2(:,1:end-1)); % dimension of feature dataset
all_mean = mean(all_feat_2(:,1:end-1));% mean of features vectors
all_no_mean = all_feat_2(:,1:end-1) - repmat(all_mean,all_sz(1),1); % features dataset normalized by mean = 0 

vat = 0.95 %variance of 95%
fs = 20 % sampling frequency 

X = all_no_mean
[l, c] = size(X) %dimention of matrix
t = [0:l-1]/fs %time
X1 = []
for i = 1:c
    X1(:,i) = X(:,i)/std(X(:,i));
end
X = X1;

S = corr(X) % correlation matrix
[e, lambda] = eig(S) % eigenvectors and eigenvalues
lambda = diag(lambda); % main diagonal of eigenvalues 

pca_perc = (lambda*100)/sum(lambda); % turning the sum of eigenvalues to 100% 
pca_n = find(cumsum(pca_perc)>= 95); % choosing the components that explain 95% of variance

%Generating principal components
cp = X * e;
pca_cp = [cp(:,1:pca_n) all_feat_2(:,28)]; 
cv = cvpartition(size(pca_cp,1),'HoldOut',0.3);
idx = cv.test;
% division of training and test data
PCA_training = pca_cp(~idx,:);
PCA_testing  = pca_cp(idx,:);
pscoreTraing95 = PCA_training(:,1:end-1)
pscoreTraing95 = PCA_testing(:,1:end-1)


%verifying if the correlation of the components is 0 (they should be)
corr_12 = corr(cp(:,1),cp(:,2))
corr_23 = corr(cp(:,2),cp(:,3))
corr_13 = corr(cp(:,1),cp(:,3))

%verifying if the dot product is 0 (they should be)
dot_12 = dot(e(:,1),e(:,2)')
dot_23 = dot(e(:,2),e(:,3)')
dot_13 = dot(e(:,1),e(:,3)')


%% Decision Trees
dt = fitctree(pscoreTrain95,trn_c);
rng(10)
dt1 = fitctree(pscoreTrain95,trn_c,'MaxNumSplits',7,'CrossVal','on','Leaveout','on','PredictorNames',{'PC1','PC2','PC3'});
view(dt1.Trained{1},'Mode','graph')
dt_rloss = dt.resubLoss;
rng(10)
dt_cv = crossval(dt,'Leaveout','on');
dt_kloss = kfoldLoss(dt_cv);
dt_accu_cv = 1-dt_kloss;
[dt_cvpred,dt_cvscores] = kfoldPredict(dt_cv);
[dtcv_CM, dtcv_ordCM] = confusionmat(dt_cvpred,trn_c);
dt_fit = predict(dt,pscoreTest95);
dt_CM = confusionmat(dt_fit,tst_c);
dt_accu = sum(diag(dt_CM))/ sum(dt_CM(:));


%% Naive Bayes
nb = fitcnb(pscoreTrain95,trn_c); %,'ClassNames',{'TP','TE','DS'}
nb_rloss = nb.resubLoss;
rng(10)
nb_cv = crossval(nb,'Leaveout','on');
nb_kloss = kfoldLoss(nb_cv);
nb_accu_cv = 1-nb_kloss;
[nb_cvpred,nb_cvscores] = kfoldPredict(nb_cv);
[nbcv_CM, nbcv_ordCM] = confusionmat(nb_cvpred,trn_c);
nb_fit = predict(nb,pscoreTest95);
nb_CM = confusionmat(nb_fit,tst_c);
nb_accu = sum(diag(nb_CM))/ sum(nb_CM(:));


%% KNN
knn = fitcknn(pscoreTrain95,trn_c,'NumNeighbors',7,'Distance','Mahalanobis');

%cvp = cvpartition(numObs,'LeaveOut')
%cvError = crossval('mcr',trn,trn_c,'Predfun',@classf,'Partition',cv_k)
%err = crossval('mse',trn,trn_c,'Predfun',@regf)
knn_rloss = knn.resubLoss;
rng(10)
knn_cv = crossval(knn,'Leaveout','on');
%knn_cv = crossval(knn,'KFold',5);
knn_kloss = kfoldLoss(knn_cv);
knn_accu_cv = 1-knn_kloss;
%matriz de confusão do cross-validated model
[knn_cvpred,knn_cvscores] = kfoldPredict(knn_cv);
[kcv_CM, kcv_ordCM] = confusionmat(knn_cvpred,trn_c);
knn_fit = predict(knn,pscoreTest95);
[knn_CM,knn_ordCM] = confusionmat(tst_c,knn_fit);
knn_accu = sum(diag(knn_CM))/ sum(knn_CM(:));

%para saber qual é o melhor número de vizinhos
knn_accu_all = [];
for i = 5:40
    knn = fitcknn(pscoreTrain95,trn_c,'NumNeighbors',i,'Distance','Mahalanobis');
    knn_fit = predict(knn,pscoreTest95);
    knn_accu_all(i) = mean(tst_c == knn_fit); 
end
figure 
plot(knn_accu_all,'-o')
xlabel('k')
ylabel('accuracy')
title('accuracy x k neighbours')


%% Random Forest

rf1 = fitensemble(pscoreTrain95,trn_c,'AdaBoostM2',100,'Tree','PredictorNames',{'PC1','PC2','PC3'})
rf2 = fitensemble(pscoreTrain95,trn_c,'AdaBoostM2',100,'Tree','CrossVal','on','Leaveout','on')
rng(10)
%{
rf1_cv = crossval(rf1,'Leaveout','on');
rf1_kloss = kfoldLoss(rf1);
fr1_accu_cv = 1-rf1_kloss;
[rf1_cvpred,rf1_cvscores] = kfoldPredict(rf1_cv);
[rf1cv_CM, rf1cv_ordCM] = confusionmat(rf1_cvpred,trn_c);
%}
rf1_kloss = kfoldLoss(rf2);
rf1_accu_cv = 1-rf1_kloss;
[rf1_cvpred,rf1_cvscores] = kfoldPredict(rf2);
[rf1cv_CM, rf1cv_ordCM] = confusionmat(rf1_cvpred,trn_c);

rf1_fit = predict(rf1,pscoreTest95);
rf1_CM = confusionmat(rf1_fit,tst_c);
rf1_accu = sum(diag(rf1_CM))/ sum(rf1_CM(:));

%% Trees of Random Forest
rng(10)
view(rf1.Trained{1},'Mode','graph')
view(rf1.Trained{6},'Mode','graph')
%view(rf1.Trained{85},'Mode','graph')

%% SVM
%MATLAB Classification learner app was used to create the model MSVM setting all the particularities automatically (gaussian kernel function and one-vs-all method)
SVM_rloss = MSVM.resubLoss;
rng(10)
SVM_cv = crossval(MSVM,'Leaveout','on');
SVM_kloss = kfoldLoss(SVM_cv);
SVM_accu_cv = 1-SVM_kloss;
[SVM_cvpred,SVM_cvscores] = kfoldPredict(SVM_cv);
[SVMcv_CM, SVMclassifiercv_ordCM] = confusionmat(SVM_cvpred,trn_c);
SVM_fit = predict(MSVM,pscoreTest95);
SVM_CM = confusionmat(SVM_fit,tst_c);
SVM_accu = sum(diag(SVM_CM))/ sum(SVM_CM(:));


%% precision, recall, F1-Measure, accuracy and confusion matrix

knn_precision = sum(diag(knn_CM))/ (sum(diag(knn_CM))+knn_CM(1,2)+knn_CM(1,3)+knn_CM(2,3));
dt_precision = sum(diag(dt_CM))/ (sum(diag(dt_CM))+dt_CM(1,2)+dt_CM(1,3)+dt_CM(2,3));
nb_precision = sum(diag(nb_CM))/ (sum(diag(nb_CM))+nb_CM(1,2)+nb_CM(1,3)+nb_CM(2,3));
rf1_precision = sum(diag(rf1_CM))/ (sum(diag(rf1_CM))+rf1_CM(1,2)+rf1_CM(1,3)+rf1_CM(2,3));
SVM_precision = sum(diag(SVM_CM))/ (sum(diag(SVM_CM))+SVM_CM(1,2)+SVM_CM(1,3)+SVM_CM(2,3));
Comparison_precision = [knn_precision dt_precision nb_precision rf1_precision SVM_precision];
Compare_Table_Precision = array2table(Comparison_precision,'VariableNames',{'KNN','DT','NB','RF','SVM'},'RowNames',{'Precision'})

knn_recall = sum(diag(knn_CM))/ (sum(diag(knn_CM))+knn_CM(2,1)+knn_CM(3,1)+knn_CM(3,2));
dt_recall = sum(diag(dt_CM))/ (sum(diag(dt_CM))+dt_CM(2,1)+dt_CM(3,1)+dt_CM(3,2));
nb_recall = sum(diag(nb_CM))/ (sum(diag(nb_CM))+nb_CM(2,1)+nb_CM(3,1)+nb_CM(3,2));
rf1_recall = sum(diag(rf1_CM))/ (sum(diag(rf1_CM))+rf1_CM(2,1)+rf1_CM(3,1)+rf1_CM(3,2));
SVM_recall = sum(diag(SVM_CM))/ (sum(diag(SVM_CM))+SVM_CM(2,1)+SVM_CM(3,1)+SVM_CM(3,2));
Comparison_recall = [knn_recall dt_recall nb_recall rf1_recall SVM_recall];
Compare_Table_Recall = array2table(Comparison_recall,'VariableNames',{'KNN','DT','NB','RF','SVM'},'RowNames',{'Recall'})

knn_F1 = 2*knn_precision*knn_recall/(knn_precision + knn_recall);
dt_F1 = 2*dt_precision*dt_recall/(dt_precision + dt_recall);
nb_F1 = 2*nb_precision*nb_recall/(nb_precision + nb_recall);
rf1_F1 = 2*rf1_precision*rf1_recall/(rf1_precision + rf1_recall);
SVM_F1 = 2*SVM_precision*SVM_recall/(SVM_precision + SVM_recall);
Comparison_F1 = [knn_F1 dt_F1 nb_F1 rf1_F1 SVM_F1];
Compare_Table_F1 = array2table(Comparison_F1,'VariableNames',{'KNN','DT','NB','RF','SVM'},'RowNames',{'F1'})


%comparing the machine learning techniques

CV_Comparison = [knn_accu_cv dt_accu_cv nb_accu_cv rf1_accu_cv SVM_accu_cv];
Comparison = [knn_accu dt_accu nb_accu rf1_accu SVM_accu];
CV_Compare_Table_Accu = array2table(CV_Comparison,'VariableNames',{'KNN_CV','DT_CV','NB_CV','RF_CV','SVM_CV'},'RowNames',{'Acuracia'})
Compare_Table_Accu = array2table(Comparison,'VariableNames',{'KNN','DT','NB','RF','SVM'},'RowNames',{'Acuracia'})

knn_CM 
dt_CM 
nb_CM 
rf1_CM
SVM_CM 





