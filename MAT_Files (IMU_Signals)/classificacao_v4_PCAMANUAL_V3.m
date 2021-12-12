clc
%{
features:
    desvio padrão no tempo
    média no tempo
    amplitude no tempo
    pico da densidade em frequencia
    frequencia dominante
    segundo pico da densidade em frequencia    
    segunda frequencia dominante
    correlação cruzada
%}

% all_feat = [astd_X; astd_Y; astd_Z; amean_X; amean_Y; amean_Z; aAmp_X; aAmp_Y; aAmp_Z; aAmp_SpX; aAmp_SpY; 
%     aAmp_SpZ; af_dom1X; af_dom1Y; af_dom1Z; aSp_2ndX; aSp_2ndY; aSp_2ndZ; af_dom2X; af_dom2Y; af_dom2Z;
%     aAmp_CRXY; aAmp_CRXZ; aAmp_CRYZ; aInd_CRXY; aInd_CRXZ; aInd_CRYZ]';

all_feat = [gstd_X; gstd_Y; gstd_Z; gmean_X; gmean_Y; gmean_Z; gAmp_X; gAmp_Y; gAmp_Z; gAmp_SpX; gAmp_SpY; 
    gAmp_SpZ; gf_dom1X; gf_dom1Y; gf_dom1Z; gSp_2ndX; gSp_2ndY; gSp_2ndZ; gf_dom2X; gf_dom2Y; gf_dom2Z;
    gAmp_CRXY; gAmp_CRXZ; gAmp_CRYZ; g_tcorr_XY; g_tcorr_XZ; g_tcorr_YZ]';

%{
all_feat_total = [gstd_X; gstd_Y; gstd_Z; gmean_X; gmean_Y; gmean_Z; gAmp_X; gAmp_Y; gAmp_Z; gAmp_SpX; gAmp_SpY; 
    gAmp_SpZ; gf_dom1X; gf_dom1Y; gf_dom1Z; gSp_2ndX; gSp_2ndY; gSp_2ndZ; gf_dom2X; gf_dom2Y; gf_dom2Z;
    gAmp_CRXY; gAmp_CRXZ; gAmp_CRYZ; g_tcorr_XY; g_tcorr_XZ; g_tcorr_YZ; astd_X;astd_Y; astd_Z; amean_X; amean_Y; amean_Z; 
    aAmp_X; aAmp_Y; aAmp_Z; aAmp_SpX; aAmp_SpY; aAmp_SpZ; af_dom1X; af_dom1Y; af_dom1Z; aSp_2ndX; aSp_2ndY; aSp_2ndZ;
    af_dom2X; af_dom2Y; af_dom2Z; aAmp_CRXY; aAmp_CRXZ; aAmp_CRYZ; a_tcorr_XY; a_tcorr_XZ; a_tcorr_YZ]';
%}
    
all_feat_X = [astd_X; amean_X; aAmp_X; aAmp_SpX; af_dom1X; aSp_2ndX; af_dom2X]';

%class = [0 0 0 0 1 1 1 1 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 2 2 0 0 2 2 0 2 2 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 2 2 0 0 0 0 0 0 1 1 0 0 2 2 2 2 0 0 0 2 2 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0];
c_012 = csvread('class_patients.csv');
%all_feat_012 = csvread('all_feat_012.csv')

Dados_HUCFF = [all_feat, c_012];
%Dados_HUCFF_Total = [all_feat_total, c_012];
csvwrite('Dados_HUCFF.csv',Dados_HUCFF)
%csvwrite('Dados_HUCFF_Completo.csv',Dados_HUCFF_Total)

%csvwrite('AF.csv',all_feat)
%figure
%plotmatrix(all_feat(1:11,:))
all_feat_table =  array2table(all_feat, 'VariableNames', {'gstd_X', 'gstd_Y', 'gstd_Z', 'gmean_X', 'gmean_Y', 'gmean_Z', 'gAmp_X', 'gAmp_Y', 'gAmp_Z', 'gAmp_SpX', 'gAmp_SpY','gAmp_SpZ', 'gf_dom1X', 'gf_dom1Y', 'gf_dom1Z', 'gSp_2ndX', 'gSp_2ndY', 'gSp_2ndZ', 'gf_dom2X', 'gf_dom2Y', 'gf_dom2Z','gAmp_CRXY', 'gAmp_CRXZ', 'gAmp_CRYZ', 'g_tcorr_XY', 'g_tcorr_XZ', 'g_tcorr_YZ'});

%balanceamento de dados
class_names = unique(c_012);
PD_cases = sum(c_012==0);
ES_cases = sum(c_012==1);
DS_cases = sum(c_012==2);
Dataset = [all_feat c_012];

 %%% balanceamento dos dados por unersampling %%%
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

%
%% Divisão em treino em teste
%aleatoridade dos dados
s_feat = size(all_feat);
numObs = length(all_feat);
rng(100)
perm = randperm(numObs);
for m = 1:s_feat(2)
    for i=1:numObs
        all_feat_2(i,m) = all_feat(perm(i),m);
    end
end

TrainPart = floor(numObs*0.7); % 70% das observações
trn = all_feat_2(1:TrainPart,1:27); % Conjunto de Treinamento
tst = all_feat_2(TrainPart+1:end,1:27); % Conjunto de Teste
trn_c = all_feat_2(1:TrainPart,28); % Classes do conjunto de treinamento
tst_c = all_feat_2(TrainPart+1:end,28); % Classes do conjunto de teste 


%% PCA train

% "all_feat_2" é a base de dados, onde as linhas são as observações e as
% colunas as características + classe de cada observação 
fs = 20;
[trn_lin trn_col] = size(trn(:,1:end)); %dimensão do base dados apenas com as características (treinamento)
trn_mean = mean(trn(:,1:end));%média dos vetores de características
trn_no_mean = trn(:,1:end) - repmat(trn_mean,trn_lin,1); %base de dados normalizada a uma média 0
t_trn = [0:trn_lin-1]/fs;
trn_norm = []
for i = 1:trn_col
    trn_norm(:,i) = trn_no_mean(:,i)/std(trn_no_mean(:,i))
end
trn_S = corr(trn_norm);
[trn_e trn_lambda] = eig(trn_S);
trn_lambda = diag(trn_lambda)
trn_perc = (trn_lambda*100)/sum(trn_lambda); %fazendo a proporção para que a soma dos autovalores dê 100% 
trn_num = find(cumsum(trn_perc)>= 95); %pegando as componentes que possuem 95% da variância 
trn_cp = trn_norm * trn_e;

trn_corr_12 = corr(trn_cp(:,1),trn_cp(:,2))
trn_corr_23 = corr(trn_cp(:,2),trn_cp(:,3))
trn_corr_13 = corr(trn_cp(:,1),trn_cp(:,3))
%% PCA test

fs = 20;
[tst_lin tst_col] = size(tst(:,1:end)); %dimensão do base dados apenas com as características (treinamento)
tst_mean = mean(tst(:,1:end));%média dos vetores de características
tst_no_mean = tst(:,1:end) - repmat(tst_mean,tst_lin,1); %base de dados normalizada a uma média 0
t_tst = [0:tst_lin-1]/fs;
tst_norm = []
for i = 1:tst_col
    tst_norm(:,i) = tst_no_mean(:,i)/std(tst_no_mean(:,i))
end
tst_S = corr(tst_norm);
[tst_e tst_lambda] = eig(tst_S);
tst_lambda = abs(diag(tst_lambda))
%tst_perc = (tst_lambda*100)/sum(tst_lambda); %fazendo a proporção para que a soma dos autovalores dê 100% 
%tst_num = find(cumsum(trn_perc)>= 95); %pegando as componentes que possuem 95% da variância 
tst_cp = tst_norm * tst_e;

tst_corr_12 = corr(tst_cp(:,1),tst_cp(:,2))
tst_corr_23 = corr(tst_cp(:,2),tst_cp(:,3))
tst_corr_13 = corr(tst_cp(:,1),tst_cp(:,3))

%%

%PCA juntando infos 


pscoreTrain95 = trn_cp(:,1:trn_num); %observações de treinamento em componentes
pscoreTest95 = tst_cp(:,1:trn_num); %observações de teste em componentes
PCA_training = [pscoreTrain95 trn_c]; %observações de treinamento + classe de cada observação
PCA_testing = [pscoreTest95 tst_c];%observações de teste + classe de cada observação

%verificando se as componentes são correlacionadas
corr_12_trn = corr(PCA_training(:,1),PCA_training(:,2))
corr_23_trn = corr(PCA_training(:,2),PCA_training(:,3))
corr_13_trn = corr(PCA_training(:,1),PCA_training(:,3))

corr_12_tst = corr(PCA_testing(:,1),PCA_testing(:,2))
corr_23_tst = corr(PCA_testing(:,2),PCA_testing(:,3))
corr_13_tst = corr(PCA_testing(:,1),PCA_testing(:,3))

%% Decision Trees
tic
rng(10)
dt = fitctree(pscoreTrain95,trn_c);
elapsed_time_dt = toc
%%
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
knn = fitcknn(pscoreTrain95,trn_c,'NumNeighbors',5,'Distance','Mahalanobis');

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


%% random forest 10 tentativa
tic% meu código MATLAB elapsed_time = toc
rng(10)
rf1 = fitensemble(pscoreTrain95,trn_c,'AdaBoostM2',100,'Tree','PredictorNames',{'PC1','PC2','PC3'});
elapsed_time_rf = toc
%%

rf2 = fitensemble(pscoreTrain95,trn_c,'AdaBoostM2',100,'Tree','CrossVal','on','Leaveout','on')
rng(10)

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

%% SVM
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


%% precision, recall, F1-Measure, acurácia, matriz de confusão

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

CV_Comparison = [knn_accu_cv dt_accu_cv nb_accu_cv rf1_accu_cv SVM_accu_cv];
Comparison = [knn_accu dt_accu nb_accu rf1_accu SVM_accu];
CV_Compare_Table_Accu = array2table(CV_Comparison,'VariableNames',{'KNN_CV','DT_CV','NB_CV','RF_CV','SVM_CV'},'RowNames',{'Acuracia'})
Compare_Table_Accu = array2table(Comparison,'VariableNames',{'KNN','DT','NB','RF','SVM'},'RowNames',{'Acuracia'})

knn_CM 
dt_CM 
nb_CM 
rf1_CM
SVM_CM 





