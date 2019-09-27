    clear;
    %datalist{1}='orlraws10P';
    %datalist{2}='pixraw10P';
    %datalist{3}='warpAR10P';
    datalist{1}='warpPIE10P';
    
    %datalist{5}='TOX-171'; 
    %datalist{6}='Carcinom';
    %datalist{7}='LUNG';
    %datalist{8}='Prostate-GE';
    %datalist{9}='GLIOMA';
    
    %datalist{10}='USPS';
    %datalist{11}='COIL20';
    
    %datalist{12}='ALLAML'; 
    %datalist{13}='CLL-SUB-111';
    %datalist{14}='GLI-85';
    %datalist{15}='SMK-CAN-187';
    %datalist{16}='isolet';    
    
num=[10:10:150]; % number of selected features
% num=[50 100];
lam_ind=-6:6; % lam_ind 的取值[-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6];
% lam_ind=[-4 4];

lam=10.^lam_ind;
p_select = 10.^lam_ind;
data_select = 1:1; % 数据集编号
% data_select = [1 2 3 4 5 6 7 8 9 10 11 12 13 15 16];
    maxRed = 1000; % maxRed
    minVal = 0;
    DictSize=1; % DictSize 
    method='MRSR';
    
% 添加当前文件夹及其子文件夹到路径 %  
currentFolder = pwd;
addpath(genpath(currentFolder));

for pi = 1:length(p_select) % pi = 1:13 (13对应于lam_ind取值的个数13个)
        for k=1:length(lam) % k = 1:13
            lambda=lam(k);
            for n=1:length(lam)
                lambda1=lam(n);
                for i=1:length(data_select)
                    fprintf('(p_value,lam_value,data)==========(%d,%d,%d)\n',pi,k,i);
                    kk=data_select(i);
                    eval(['load ' datalist{kk}])
                    fprintf(datalist{kk});
                    fprintf('\n');
                    X = NormalizeFea(X);
                    sam_num=size(X,1); % sam_num: 样本数
                    feature_num=size(X,2); % feature_num：特征数

                    [indx stump]=evalute(1,X,Y,method,feature_num,lambda,lambda1,p_select(pi),DictSize);  % 得到特征得分从高到低排序好的矩阵

                    for j=1:length(num) % 取num=[10:10:150]个特征评价结果
                        feature_num=num(j); %   j=1:15对应于的 feature_num 为 10*j        
                        [rec_acc_fs(i,j),rec_clu(i,j),rebunduncy(i,j),rec_acc_clu(i,j)]=evalute_num(X,Y,feature_num,indx); % (i,j):(dataset_id,feature_Num_id)
                    % classification，NMI，       rebunduncy，    clustering
                    end  
                end
                recacc{k}=rec_acc_fs'; % recacc:一维cell   其中recacc的第k个cell记录相应的结果值，k对应lamda的取值个数。 
                recclu{k}=rec_clu'; 
                rebund{k}=rebunduncy'; 
                recacclu{k}=rec_acc_clu'; 

                pLamDataFeanumAcc{pi,k}=rec_acc_fs;  % pLamDataFeanumAcc: pi*k维cell，(pi,k)中记录 第pi个lamda，第k个lamda对应的结果值（13*13）。rec_acc_fs(i,j)中记录第i个数据集，第j个特征数取值对应的结果值。
                pLamDataFeanumClu{pi,k}=rec_clu;
                pLamDataFeanumUnd{pi,k}=rebunduncy;
                pLamDataFeanumAccclu{pi,k}=rec_acc_clu;
            end
        end
end    
for k=1:length(data_select) % k = 1:16 ; 对应于16个数据集
    s=0;
    result1=[];
    result2=[];
    result3=[];
    result4=[];
 for m=1:length(p_select)  % m = 1:13 ; 对应于lamda的13个取值 
     for n=1:length(lam) % n = 1:13 ; 
            s=s+1; % s的取值：1…n,n+1…m*n（13*13=169）；一共k轮
            result1(s)=mean(pLamDataFeanumAcc{m,n}(k,:)); % 对于pLamDataFeanumAcc的 特征数取值选项 取均值。result1(s)行向量。
            result2(s)=mean(pLamDataFeanumClu{m,n}(k,:));
            result3(s)=mean(pLamDataFeanumUnd{m,n}(k,:));
            result4(s)=mean(pLamDataFeanumAccclu{m,n}(k,:));
     end
 end  

value1(k)=max(result1); % 取出对应于第k个数据集的 result1中的最大值
value2(k)=max(result2);
value3(k)=min(result3);
value4(k)=max(result4);

% deviation1(k)=std(result1);
% deviation2(k)=std(result2);
% deviation3(k)=std(result3);
% deviation4(k)=std(result4);

end

result=[value1',value2',value3',value4']*100
% result_deviation=[deviation1',deviation2',deviation3',deviation4']*100




