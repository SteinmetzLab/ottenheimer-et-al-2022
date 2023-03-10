%%  get raw spikes at different time bins for decoding
dwindow=[-0.5 2.5];
dbinsize=0.25;
dbinEdges=dwindow(1):dbinsize:dwindow(2);
dbinTimes=dwindow(1)+dbinsize/2:dbinsize:dwindow(2)-dbinsize/2;

preRun = true;

if preRun
    load('dactivity.mat');
else
    dactivity=cell(totalNeurons,6);
    NN=0;
    for incses=1:sum(includedSessions)
        session=sessInd(incses);
        %get events
        sessionFolder=fullfile(direc,sessionSubject{session},sessionDate{session});
        cueTimes = readNPY(fullfile(sessionFolder,'cues.times.npy'));
        rewprob = readNPY(fullfile(sessionFolder,'cues.rewardProbability.npy'));
        odorset = readNPY(fullfile(sessionFolder,'cues.odorSet.npy'));
        csp1 = cueTimes(rewprob==1&odorset==1);
        csp2 = cueTimes(rewprob==1&odorset==2);
        csf1 = cueTimes(rewprob==0.5&odorset==1);
        csf2 = cueTimes(rewprob==0.5&odorset==2);
        csm1 = cueTimes(rewprob==0&odorset==1);
        csm2 = cueTimes(rewprob==0&odorset==2);
        trialNo=1:length(cueTimes);
        
        folders=dir(sessionFolder);
        probeFolderIDs = find(arrayfun(@(x) strcmp('p',folders(x).name(1)),1:length(folders)));
        if isreal(probeFolderIDs)
            for probe=1:length(probeFolderIDs)
                probeFolder=fullfile(sessionFolder,folders(probeFolderIDs(probe)).name);
                shankDir=dir(probeFolder);
                entryName={};
                entryNameWithShank={};
                for entry=1:length(shankDir)
                    entryName{entry,1}=shankDir(entry).name(1:min([length(shankDir(entry).name) 5]));
                    entryNameWithShank{entry,1}=shankDir(entry).name(1:min([length(shankDir(entry).name) 6]));
                end
                shanks=sum(strcmp(entryName,'shank'));
                shankList=[];
                for entry=1:length(shankDir)
                    if length(entryNameWithShank{entry})>5
                        shankList=cat(1,shankList,str2num(entryNameWithShank{entry}(6)));
                    end
                end
                for shank=1:shanks
                    
                    shankFolder=fullfile(probeFolder,append('shank',num2str(shankList(shank))));
                    clusterNames=readNPY(fullfile(shankFolder,'clusters.names.npy'));
                    incClusters=clusterNames(readNPY(fullfile(shankFolder,'clusters.passedQualityControl.npy')));
                    channelLocations = textscan(fopen(fullfile(shankFolder,'channels.CCFregions.txt')),'%s');
                    fclose('all');
                    channelLocations = channelLocations{1};
                    spikeTimes = readNPY(fullfile(shankFolder,'spikes.times.npy'));
                    spikeClusters = readNPY(fullfile(shankFolder,'spikes.clusterNames.npy'));
                    clusterChannel = readNPY(fullfile(shankFolder,'clusters.channels.npy'));
                    
                    %get clusters
                    for clus=1:length(incClusters)
                        
                        cluster = incClusters(clus);
                        regionName = channelLocations{clusterChannel(clusterNames==cluster)};
                        charOnly = regexp(regionName,'[c-zA-Z]','match');
                        noNums=cat(2,charOnly{:});
                        if sum(strcmp(allAreas,noNums(1:min([3 length(noNums)]))))==1
                            NN=NN+1;
                            clusterTimes = double(spikeTimes(spikeClusters==cluster))/30000;
                            
                            
                            binnedActivity=NaN(length(cueTimes),length(dbinEdges)-1);
                            numBins=length(dbinTimes);
                            for trial = 1:length(cueTimes)
                                clusterTimesRel = clusterTimes - cueTimes(trial);
                                binnedActivity(trial,:)=histcounts(clusterTimesRel,dbinEdges);
                            end
                            binnedActivity(isnan(binnedActivity))=0;binnedActivity(binnedActivity==inf)=0; %this only happens if bins go beyond end of session, so extremely rare
                            
                            %3 conditions
                            selections={};
                            selections{1,1}=csp1; %cs+
                            selections{2,1}=csp2; %cs+
                            selections{3,1}=csf1; %cs50
                            selections{4,1}=csf2; %cs50
                            selections{5,1}=csm1; %cs-
                            selections{6,1}=csm2; %cs-
                            for condition=1:length(selections)
                                sel = ismember(round(cueTimes),round(selections{condition}));
                                dactivity{NN,condition}=binnedActivity(sel,:);
                            end
                        end
                    end
                end
            end
        end
        
        
        fprintf('Session %d \n',session);
    end
    
    save('dactivity.mat','dactivity','-v7.3');
    
end

    
%% single unit decoding of cue identity
preRun=true;

if preRun==false
dts=50;
folds=5;
trains=ones(dts*6,folds);
for t=1:folds
    trains(t:5:end,t)=0;
end

odorLabels=[ones(dts,1);2*ones(dts,1);3*ones(dts,1);4*ones(dts,1);5*ones(dts,1);6*ones(dts,1)];
valueLabels=[ones(dts,1);ones(dts,1);0.5*ones(dts,1);0.5*ones(dts,1);zeros(dts,1);zeros(dts,1)];
setLabels=[ones(dts,1);2*ones(dts,1);ones(dts,1);2*ones(dts,1);ones(dts,1);2*ones(dts,1)];

tps=size(dactivity{1,1},2);
odorPerf=NaN(totalNeurons,tps);
valuePerf=NaN(totalNeurons,tps);
setPerf=NaN(totalNeurons,tps);
tic
for NN=1:totalNeurons
    
    activity=[];
    for cue=1:6
        trls=randperm(length(dactivity{NN,cue}),dts);
        activity=cat(1,activity,dactivity{NN,cue}(trls,:));
    end
    
    parfor t=1:tps
        predictionO=NaN(dts*6,1);
        predictionV=NaN(dts*6,1);
        predictionS=NaN(dts*6,1);
        for fold=1:folds
            sel=trains(:,fold)==1;
            aoi=activity(sel,t);
            if sum(aoi)>0 %only run if there are spikes
            mdl = fitcdiscr(aoi,odorLabels(sel));
            predictionO(sel==0)=predict(mdl,activity(sel==0,t));
            
            mdl = fitcdiscr(aoi,valueLabels(sel));
            predictionV(sel==0)=predict(mdl,activity(sel==0,t));            

            mdl = fitcdiscr(aoi,setLabels(sel));
            predictionS(sel==0)=predict(mdl,activity(sel==0,t));            
            end
        end
        difference=predictionO-odorLabels;
        odorPerf(NN,t)=sum(difference==0)/length(difference);
        
        difference=predictionV-valueLabels;
        valuePerf(NN,t)=sum(difference==0)/length(difference);        
        
        difference=predictionS-setLabels;
        setPerf(NN,t)=sum(difference==0)/length(difference);        
    end
    if rem(NN,10)==0 fprintf('Neuron %d \n',NN); end
end
toc
save('dactivity.mat','dactivity','odorPerf','valuePerf','setPerf','-v7.3');
end
%% plot single neuron decoding of cue identity
cats={1,2,3};
catcolors={[0 0.7 1],[0.7 0 1],[0.6 0.6 0.6],[0.05 0.25 0.45]};
vdata=[];
vcat=[];
vsess=[];
vtp=[];
odata=[];
for c=1:length(cats)
    sel=ismember(category,cats{c});
    
    subplot(1,3,1)
    hold on
    up=nanmean(valuePerf(sel,:),1)+nanste(valuePerf(sel,:),1);
    down=nanmean(valuePerf(sel,:),1)-nanste(valuePerf(sel,:),1);
    plot(dbinTimes,nanmean(valuePerf(sel,:),1),'color',catcolors{c},'linewidth',1);  
    patch([dbinTimes,dbinTimes(end:-1:1)],[up,down(end:-1:1)],catcolors{c},'EdgeColor','none');alpha(0.5);    
    plot(dwindow,[1/3 1/3],':','color','k','linewidth',0.75)
    
    title('Single unit decoding of value');
    xlim(dwindow);
    ylim([0.3 0.5]);
    
    ylabel('Decoding accuracy');
    xlabel('Time from odor onset (s)');

    vdata=cat(1,vdata,valuePerf(sel,:));
    vcat=cat(1,vcat,c*ones(sum(sel),1));
    vsess=cat(1,vsess,neuronSession(sel));
    vtp=cat(1,vtp,repmat(1:size(valuePerf,2),sum(sel),1));
    
    subplot(1,3,2)
    hold on
    up=nanmean(odorPerf(sel,:),1)+nanste(odorPerf(sel,:),1);
    down=nanmean(odorPerf(sel,:),1)-nanste(odorPerf(sel,:),1);
    plot(dbinTimes,nanmean(odorPerf(sel,:),1),'color',catcolors{c},'linewidth',1);  
    patch([dbinTimes,dbinTimes(end:-1:1)],[up,down(end:-1:1)],catcolors{c},'EdgeColor','none');alpha(0.5);    
    plot(dwindow,[1/6 1/6],':','color','k','linewidth',0.75)
    
    title('Single unit decoding of odor');
    xlim(dwindow);
    ylim([0.15 0.27]); 

    odata=cat(1,odata,odorPerf(sel,:));  

end

%value anova
vcat=repmat(vcat,1,size(valuePerf,2));
vsess=repmat(vsess,1,size(valuePerf,2));

%odor anova
[~,tbl,stats]=anovan(odata(:),{categorical(vcat(:)),categorical(vtp(:)),categorical(vsess(:))},'varnames',{'cat','time','sess'},'random',[3],'display','off',...
    'model',[1 0 0;0 1 0;0 0 1;1 1 0]);
c=multcompare(stats,'dimension',[1 2],'ctype','bonferroni');
%% population decoding of cue identity using different cue types
numN=[1 5 10 25 50 75 100 200];
pwindow=[1 2.5];
tpsel=dbinTimes>pwindow(1) & dbinTimes<pwindow(2);
reps=1000;
preRun=false;

if preRun==false

dts=50;
gamma=1;

folds=5;
trains=ones(dts*6,folds);
for t=1:folds
    trains(t:5:end,t)=0;
end

odorLabels=[ones(dts,1);2*ones(dts,1);3*ones(dts,1);4*ones(dts,1);5*ones(dts,1);6*ones(dts,1)];
valueLabels=[ones(dts,1);ones(dts,1);0.5*ones(dts,1);0.5*ones(dts,1);zeros(dts,1);zeros(dts,1)];

cats={1,2,3};
odorPerfP=NaN(reps,length(numN),length(cats));
valuePerfP=NaN(reps,length(numN),length(cats));

for c=1:length(cats)
    sel=ismember(category,cats{c});
    dactivityInc=dactivity(sel,:);
    for nn=1:length(numN)
        n2u=numN(nn);
        parfor r=1:reps
            nis=randperm(size(dactivityInc,1),n2u);
            activity=NaN(dts*6,n2u);
            for n=1:n2u
                ni=nis(n);
                nactivity=[];
                for cue=1:6
                    trls=randperm(length(dactivityInc{ni,cue}),dts);
                    nactivity=cat(1,nactivity,sum(dactivityInc{ni,cue}(trls,tpsel),2));
                end
                if max(nactivity,[],'all')>0
                nactivity=nactivity/max(nactivity,[],'all');
                end
                activity(:,n)=nactivity;
            end
            
            predictionO=NaN(dts*6,1);
            predictionV=NaN(dts*6,1);
            for fold=1:folds
                sel=trains(:,fold)==1;
                aoi=activity(sel,:);
                aot=activity(sel==0,:);
                if any(sum(aoi,1)==0) %only include neuron if there are spikes
                    aot=activity(sel==0,sum(aoi,1)>0);
                    aoi=activity(sel,sum(aoi,1)>0);
                end
                
                mdl = fitcdiscr(aoi,odorLabels(sel),'gamma',gamma);
                predictionO(sel==0)=predict(mdl,aot);
                
                mdl = fitcdiscr(aoi,valueLabels(sel),'gamma',gamma);
                predictionV(sel==0)=predict(mdl,aot);
                
            end
            difference=predictionO-odorLabels;
            odorPerfP(r,nn,c)=sum(difference==0)/length(difference);
            
            difference=predictionV-valueLabels;
            valuePerfP(r,nn,c)=sum(difference==0)/length(difference);
            
    
        end
        fprintf('Category %d, %d neurons \n',c,n2u);
    end
end

save('dactivity.mat','dactivity','odorPerf','valuePerf','setPerf','odorPerfP','valuePerfP','-v7.3');

end

%% plot population decoding of cue identity (by neuron category)
figure;
cats={1,2,3};
catcolors={[0 0.7 1],[0.7 0 1],[0.6 0.6 0.6],[0.05 0.25 0.45]};

subplot(1,3,1)
hold on
for c=1:length(cats)
    errorbar([1:length(numN)]-0.3+0.15*c,nanmean(squeeze(valuePerfP(:,:,c)),1),nanstd(squeeze(valuePerfP(:,:,c)),1),'o','color',catcolors{c},'linewidth',1);
end
plot([0.5 length(numN)+0.5],[1/3 1/3],':','color','k','linewidth',0.75)
ylabel('Decoding accuracy');
xlabel('Neurons used');
ylim([0 1.1]);
xlim([0.5 length(numN)+0.5]);
yticks(0:0.25:1);
xticks(1:length(numN));
xticklabels(numN);
title('Value');

%pvals
pvals={};
pvalCh=NaN(length(cats)*2,length(numN));
for n=1:length(numN)
    for c1=1:length(cats)
        for c2=1:length(cats)
            pvals{1,n}(c1,c2)=(sum(squeeze(valuePerfP(:,n,c1))<=squeeze(valuePerfP(:,n,c2))','all')+1)/(size(valuePerfP,1)^2+1);
        end
        pvalCh(c1,n)=(sum(squeeze(valuePerfP(:,n,c1))<=(1/3))+1)/(size(valuePerfP,1)+1);
        
    end
end

subplot(1,3,2)
hold on
for c=1:length(cats)
    errorbar([1:length(numN)]-0.25+0.1*c,nanmean(squeeze(odorPerfP(:,:,c)),1),nanstd(squeeze(valuePerfP(:,:,c)),1),'o','color',catcolors{c},'linewidth',1);
end
plot([0.5 length(numN)+0.5],[1/6 1/6],':','color','k','linewidth',0.75)
ylim([0 1.1]);
xlim([0.5 length(numN)+0.5]);
yticks(0:0.25:1);
xticks(1:length(numN));
xticklabels(numN);
title('Odor');

%pvals
for n=1:length(numN)
    for c1=1:length(cats)
        for c2=1:length(cats)
            pvals{2,n}(c1,c2)=(sum(squeeze(odorPerfP(:,n,c1))<=squeeze(odorPerfP(:,n,c2))','all')+1)/(size(valuePerfP,1)^2+1);
        end
        pvalCh(c1+3,n)=(sum(squeeze(odorPerfP(:,n,c1))<=(1/6))+1)/(size(valuePerfP,1)+1);

    end
end


%% population decoding of value using different cue types
numN=[1 5 10 25 50 75 100 200];
pwindow=[1 2.5];
tpsel=dbinTimes>pwindow(1) & dbinTimes<pwindow(2);
reps=1000;
preRun=false;

if preRun==false

dts=50;
gamma=0.01;

opts.alpha=0.5;
glmopts=glmnetSet(opts);

folds=6;
tests=zeros(dts*6,folds);
for f=1:folds
tests((f-1)*dts+1:f*dts,f)=1;
end

trains1=zeros(dts*6,1);trains1(51:100)=1;trains1(126:175)=1;trains1(226:275)=1;
trains2=zeros(dts*6,1);trains2(1:50)=1;trains2(126:175)=1;trains2(226:275)=1;
trains3=zeros(dts*6,1);trains3(26:75)=1;trains3(151:200)=1;trains3(226:275)=1;
trains4=zeros(dts*6,1);trains4(26:75)=1;trains4(101:150)=1;trains4(226:275)=1;
trains5=zeros(dts*6,1);trains5(26:75)=1;trains5(126:175)=1;trains5(251:300)=1;
trains6=zeros(dts*6,1);trains6(26:75)=1;trains6(126:175)=1;trains6(201:250)=1;
trains=[trains1 trains2 trains3 trains4 trains5 trains6];

valueLabels=[ones(dts,1);ones(dts,1);0.5*ones(dts,1);0.5*ones(dts,1);zeros(dts,1);zeros(dts,1)];

cats={1,2,3};
valuePerfP=NaN(reps,length(numN),length(cats));
valueCueP=NaN(reps,length(numN),length(cats),3);
for c=1:length(cats)
    sel=ismember(category,cats{c});
    dactivityInc=dactivity(sel,:);
    for nn=1:length(numN)
        n2u=numN(nn);
        parfor r=1:reps
            nis=randperm(size(dactivityInc,1),n2u);
            activity=NaN(dts*6,n2u);
            for n=1:n2u
                ni=nis(n);
                nactivity=[];
                for cue=1:6
                    trls=randperm(length(dactivityInc{ni,cue}),dts);
                    nactivity=cat(1,nactivity,sum(dactivityInc{ni,cue}(trls,tpsel),2));
                end
                if max(nactivity,[],'all')>0
                end
                activity(:,n)=nactivity;
            end
            
            predictionV=NaN(dts*6,1);
            for fold=1:folds
                trainsel=trains(:,fold)==1;
                testsel=tests(:,fold)==1;
                aoi=activity(trainsel,:);
                aot=activity(testsel,:);
                if any(sum(aoi,1)==0) %only include neuron if there are spikes
                    aot=activity(testsel,sum(aoi,1)>0);
                    aoi=activity(trainsel,sum(aoi,1)>0);
                end
                
                fit = glmnet(aoi,valueLabels(trainsel),'gaussian',glmopts);
                prediction = glmnetPredict(fit,aot,[],'response');
                predictionV(testsel)=prediction(:,end);

                
            end
            
            catPrediction=NaN(length(predictionV),1);
            catPrediction(predictionV<0.25 & predictionV>=-0.25)=0;
            catPrediction(predictionV<0.75 & predictionV>=0.25)=0.5;
            catPrediction(predictionV<1.25 & predictionV>=0.75)=1;
            
            difference=catPrediction-valueLabels;
            valuePerfP(r,nn,c)=sum(difference==0)/length(difference);
            
            for cue=1:3
                valueCueP(r,nn,c,cue)=mean(predictionV((cue-1)*100+1:cue*100));
            end
    
        end
        fprintf('Category %d, %d neurons \n',c,n2u);
    end
end

save('dactivityV.mat','dactivity','valuePerf','valuePerfP','-v7.3');

end

%% plot population decoding of value (by neuron category)
figure;
cats={1,2,3};
catcolors={[0 0.7 1],[0.7 0 1],[0.6 0.6 0.6],[0.05 0.25 0.45]};

subplot(1,3,1)
hold on
for c=1:length(cats)
    errorbar([1:length(numN)]-0.3+0.15*c,nanmean(squeeze(valuePerfP(:,:,c)),1),nanstd(squeeze(valuePerfP(:,:,c)),1),'o','color',catcolors{c},'linewidth',1);
end
plot([0.5 length(numN)+0.5],[1/3 1/3],':','color','k','linewidth',0.75)
ylabel('Decoding accuracy');
xlabel('Neurons used');
ylim([0 1.1]);
xlim([0.5 length(numN)+0.5]);
yticks(-1:0.5:1);
xticks(1:length(numN));
xticklabels(numN);
title('Value');

   subplot(2,3,2)
    hold on
for c=1:length(cats)
    errorbar([1:length(numN)]-0.3+0.15*c,nanmean(squeeze(valueCueP(:,:,c,1)),1),nanstd(squeeze(valueCueP(:,:,c,1)),1),'o','color',catcolors{c},'linewidth',1);
end
plot([0.5 length(numN)+0.5],[1 1],':','color',colors{1},'linewidth',0.75)
plot([0.5 length(numN)+0.5],[0.5 0.5],':','color',colors{2},'linewidth',0.75)

ylabel('Predicted cue value');
ylim([0.4 1.1]);
xlim([0.5 length(numN)+0.5]);
yticks(-1:0.5:1);
xticks([]);
xticklabels(numN);
title('CS+');

   subplot(2,3,5)
    hold on
for c=1:length(cats)
    errorbar([1:length(numN)]-0.3+0.15*c,nanmean(squeeze(valueCueP(:,:,c,3)),1),nanstd(squeeze(valueCueP(:,:,c,3)),1),'o','color',catcolors{c},'linewidth',1);
end
plot([0.5 length(numN)+0.5],[0 0],':','color',colors{3},'linewidth',0.75)
plot([0.5 length(numN)+0.5],[0.5 0.5],':','color',colors{2},'linewidth',0.75)

xlabel('Neurons used');
ylim([-0.1 0.6]);
xlim([0.5 length(numN)+0.5]);
yticks(-1:0.5:1);
xticks(1:length(numN));
xticklabels(numN);
title('CS-');

%pvals
pvals={};
pvalCh=NaN(length(cats)*2,length(numN));
for n=1:length(numN)
    for c1=1:length(cats)
        for c2=1:length(cats)
            pvals{1,n}(c1,c2)=(sum(squeeze(valuePerfP(:,n,c1))<=squeeze(valuePerfP(:,n,c2))','all')+1)/(size(valuePerfP,1)^2+1);
        end
        pvalCh(c1,n)=(sum(squeeze(valuePerfP(:,n,c1))<=(1/3))+1)/(size(valuePerfP,1)+1);
        
    end
end

%% population decoding of value using different regions
n2u=5;
pwindows={[-0.5 0],[1 2.5]};
reps=1000;
preRun=false;

regions={'ALM','ACA','FRP','PL','ILA','ORB','DP','TTd','AON'};
neuronRegionOlf=neuronRegionAdj;
neuronRegionOlf(ismember(neuronRegionAdj,{'EPd','PIR'}))={'OLF'};

if preRun==false


folds=6;
tests=zeros(dts*6,folds);
for f=1:folds
tests((f-1)*dts+1:f*dts,f)=1;
end

trains1=zeros(dts*6,1);trains1(51:100)=1;trains1(126:175)=1;trains1(226:275)=1;
trains2=zeros(dts*6,1);trains2(1:50)=1;trains2(126:175)=1;trains2(226:275)=1;
trains3=zeros(dts*6,1);trains3(26:75)=1;trains3(151:200)=1;trains3(226:275)=1;
trains4=zeros(dts*6,1);trains4(26:75)=1;trains4(101:150)=1;trains4(226:275)=1;
trains5=zeros(dts*6,1);trains5(26:75)=1;trains5(126:175)=1;trains5(251:300)=1;
trains6=zeros(dts*6,1);trains6(26:75)=1;trains6(126:175)=1;trains6(201:250)=1;
trains=[trains1 trains2 trains3 trains4 trains5 trains6];

valueLabels=[ones(dts,1);ones(dts,1);0.5*ones(dts,1);0.5*ones(dts,1);zeros(dts,1);zeros(dts,1)];

tps=size(dactivity{1,1},2);

valuePerfR=NaN(reps,length(regions),length(pwindows));

for tp=1:length(pwindows)
    tpsel=dbinTimes>pwindows{tp}(1) & dbinTimes<pwindows{tp}(2);
    for c=1:length(regions)
        sel=ismember(category,[1]) & ismember(neuronRegionOlf,regions(c));
        dactivityInc=dactivity(sel,:);
        
        parfor r=1:reps
            nis=randsample(size(dactivityInc,1),n2u,'true');
            activity=NaN(dts*6,n2u);
            for n=1:n2u
                ni=nis(n);
                nactivity=[];
                for cue=1:6
                    trls=randperm(length(dactivityInc{ni,cue}),dts);
                    nactivity=cat(1,nactivity,sum(dactivityInc{ni,cue}(trls,tpsel),2));
                end
                if max(nactivity,[],'all')>0
                    nactivity=nactivity/max(nactivity,[],'all');
                end
                activity(:,n)=nactivity;
            end
            
            predictionV=NaN(dts*6,1);
            for fold=1:folds
                trainsel=trains(:,fold)==1;
                testsel=tests(:,fold)==1;
                aoi=activity(trainsel,:);
                aot=activity(testsel,:);
                if any(sum(aoi,1)==0) %only include neuron if there are spikes
                    aot=activity(testsel,sum(aoi,1)>0);
                    aoi=activity(trainsel,sum(aoi,1)>0);
                end

                fit = glmnet(aoi,valueLabels(trainsel),'gaussian',glmopts);
                prediction = glmnetPredict(fit,aot,[],'response');
                predictionV(testsel)=prediction(:,end);                
                
            end
            
            catPrediction=NaN(length(predictionV),1);
            catPrediction(predictionV<0.25 & predictionV>=-0.25)=0;
            catPrediction(predictionV<0.75 & predictionV>=0.25)=0.5;
            catPrediction(predictionV<1.25 & predictionV>=0.75)=1;
            
            difference=catPrediction-valueLabels;
            valuePerfR(r,c,tp)=sum(difference==0)/length(difference);
            
        end
        fprintf('TP %d, region %d \n',tp,c);
    end
end
save('dactivityV.mat','dactivity','valuePerf','valuePerfP','valuePerfR','-v7.3');

end
%% plot population decoding of value (by region)
figure;
catcolors={[0 0 0],[0.35 0.35 1]};
catcolors={[0 0 0],[0 0.7 1]};

subplot(1,3,1)
hold on
for tp=1:length(pwindows)
errorbar(1:length(regions),nanmean(squeeze(valuePerfR(:,:,tp)),1),nanstd(squeeze(valuePerfR(:,:,tp)),1),'o','color',catcolors{tp},'linewidth',1);
end
plot([0.5 length(regions)+0.5],[1/3 1/3],':','color','k','linewidth',0.75)
ylabel('Decoding accuracy');
xlabel('Regions');
ylim([0 1.1]);
xlim([0.5 length(regions)+0.5]);
yticks(0:0.25:1);
xticks(1:length(regions));
xticklabels(regions);
xtickangle(45);
title('Value');

%pvals
pvalsV=[];
pvalB=NaN(2,length(regions));
for c1=1:length(regions)
    for c2=1:length(regions)
        pvalsV(c1,c2)=(sum(squeeze(valuePerfR(:,c1,2))<=squeeze(valuePerfR(:,c2,2))','all')+1)/(size(valuePerfR,1)^2+1);
    end
    pvalsB(1,c1)=(sum(squeeze(valuePerfR(:,c1,2))<=squeeze(valuePerfR(:,c1,1))','all')+1)/(size(valuePerfR,1)^2+1);
    
    
end


%% population decoding of value using different region groups
n2u=25;
pwindows={[-0.5 0],[1 2.5]};
reps=1000;
preRun=false;

reggrps=[1:3];
if preRun==false

dts=50;
gamma=1;

folds=6;
tests=zeros(dts*6,folds);
for f=1:folds
tests((f-1)*dts+1:f*dts,f)=1;
end

trains1=zeros(dts*6,1);trains1(51:100)=1;trains1(126:175)=1;trains1(226:275)=1;
trains2=zeros(dts*6,1);trains2(1:50)=1;trains2(126:175)=1;trains2(226:275)=1;
trains3=zeros(dts*6,1);trains3(26:75)=1;trains3(151:200)=1;trains3(226:275)=1;
trains4=zeros(dts*6,1);trains4(26:75)=1;trains4(101:150)=1;trains4(226:275)=1;
trains5=zeros(dts*6,1);trains5(26:75)=1;trains5(126:175)=1;trains5(251:300)=1;
trains6=zeros(dts*6,1);trains6(26:75)=1;trains6(126:175)=1;trains6(201:250)=1;
trains=[trains1 trains2 trains3 trains4 trains5 trains6];

valueLabels=[ones(dts,1);ones(dts,1);0.5*ones(dts,1);0.5*ones(dts,1);zeros(dts,1);zeros(dts,1)];

valuePerfG=NaN(reps,length(reggrps),length(pwindows));

for tp=1:length(pwindows)
    tpsel=dbinTimes>pwindows{tp}(1) & dbinTimes<pwindows{tp}(2);
    for c=1:length(reggrps)
        sel=ismember(category,[1]) & ismember(neuronGroup,reggrps(c));
        dactivityInc=dactivity(sel,:);
        
        parfor r=1:reps
            nis=randsample(size(dactivityInc,1),n2u,'true');
            activity=NaN(dts*6,n2u);
            for n=1:n2u
                ni=nis(n);
                nactivity=[];
                for cue=1:6
                    trls=randperm(length(dactivityInc{ni,cue}),dts);
                    nactivity=cat(1,nactivity,sum(dactivityInc{ni,cue}(trls,tpsel),2));
                end
                if max(nactivity,[],'all')>0
                    nactivity=nactivity/max(nactivity,[],'all');
                end
                activity(:,n)=nactivity;
            end
            
            predictionV=NaN(dts*6,1);
            for fold=1:folds
                trainsel=trains(:,fold)==1;
                testsel=tests(:,fold)==1;
                aoi=activity(trainsel,:);
                aot=activity(testsel,:);
                if any(sum(aoi,1)==0) %only include neuron if there are spikes
                    aot=activity(testsel,sum(aoi,1)>0);
                    aoi=activity(trainsel,sum(aoi,1)>0);
                end
                
                fit = glmnet(aoi,valueLabels(trainsel),'gaussian',glmopts);
                prediction = glmnetPredict(fit,aot,[],'response');
                predictionV(testsel)=prediction(:,end);                
                
            end
            
            catPrediction=NaN(length(predictionV),1);
            catPrediction(predictionV<0.25 & predictionV>=-0.25)=0;
            catPrediction(predictionV<0.75 & predictionV>=0.25)=0.5;
            catPrediction(predictionV<1.25 & predictionV>=0.75)=1;
            
            difference=catPrediction-valueLabels;
            valuePerfG(r,c,tp)=sum(difference==0)/length(difference);
            
        end
        fprintf('TP %d, region %d \n',tp,c);
    end
end
save('dactivityV.mat','dactivity','valuePerf','valuePerfP','valuePerfR','valuePerfG','-v7.3');

end
%% plot population decoding of value (by region group)
figure;
catcolors={[0 0 0],[0.35 0.35 1]};
catcolors={[0 0 0],[0 0.7 1]};

subplot(1,6,1)
hold on
for tp=1:length(pwindows)
errorbar(1:length(reggrps),nanmean(squeeze(valuePerfG(:,:,tp)),1),nanstd(squeeze(valuePerfG(:,:,tp)),1),'o','color',catcolors{tp},'linewidth',1);
end
plot([0.5 length(reggrps)+0.5],[1/3 1/3],':','color','k','linewidth',0.75)
ylabel('Decoding accuracy');
xlabel('Region groups');
ylim([0 1.1]);
xlim([0.5 length(reggrps)+0.5]);
yticks(0:0.25:1);
xticks(1:length(reggrps));
xticklabels(regionGroupNames(reggrps));
xtickangle(45);
title('Value');

%pvals
pvalsV=[];
for c1=1:length(reggrps)
    for c2=1:length(reggrps)
        pvalsV(c1,c2)=(sum(squeeze(valuePerfG(:,c1,2))<=squeeze(valuePerfG(:,c2,2))','all')+1)/(size(valuePerfG,1)^2+1);
    end
end

