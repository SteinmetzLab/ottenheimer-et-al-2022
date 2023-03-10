% Script to analyze electrophysiology data for value coding
githubDir = 'D:\GitHub';
addpath(genpath(fullfile(githubDir, 'npy-matlab'))) %to load .npy files
addpath(genpath(fullfile(githubDir, 'steinmetz-et-al-2019'))) %for reduced rank
direc = 'D:\GitHub\ottenheimer-et-al-2022\Neuropixels'; %where neuropixels data exists
addpath(direc);
%% find all sessions

subjectFolders = dir(direc);
nameLengths = arrayfun(@(x) length(subjectFolders(x).name),1:length(subjectFolders));
subjectFolders = subjectFolders(nameLengths>5);

sessionSubject={};
sessionDate={};
sessionNumber=0;
for sn=1:length(subjectFolders)
    dateFolders = dir(fullfile(direc,subjectFolders(sn).name));
    nameLengths = arrayfun(@(x) length(dateFolders(x).name),1:length(dateFolders));
    dateFolders = dateFolders(nameLengths>5);
    for d=1:length(dateFolders)
        sessionNumber=sessionNumber+1;
        sessionSubject{sessionNumber,1}=subjectFolders(sn).name;
        sessionDate{sessionNumber,1}=dateFolders(d).name;
    end
end

%% count how many neurons

%region categories
tic
NN=0;
NP=0;
probeNum=0;
probeNeurons=[];
neuronRegion={};
neuronRegionAdj={};
neuronSubregion={};
neuronGroup=[];
neuronSubject={};
neuronDate={};
neuronCluster=[];
neuronSession=[];
neuronProbe=[];
neuronRegionID=[];
neuronHz=[];
neuronXYZ=[];
includedSessions=false(sessionNumber,1);
sessionNeurons=NaN(sessionNumber,1);


motorCortex = {'MOs'};
olfactoryAreas = {'AON','DP','EPd','OLF','TTd','PIR'};
assocCortex = {'ACA','ILA','ORB','PL','FRP'};
striatalAreas = {'ACB','CP'};
allAreas = cat(2,motorCortex,olfactoryAreas,assocCortex,striatalAreas);
regionGroupNames={'Motor','PFC','Olfactory','Striatum'};
regionGroups={motorCortex;assocCortex;olfactoryAreas;striatalAreas};

tic
for session=1:sessionNumber
    NNstart=NN;
    sessionFolder=fullfile(direc,sessionSubject{session},sessionDate{session});
    folders=dir(sessionFolder);
    probeFolderIDs = find(arrayfun(@(x) strcmp('p',folders(x).name(1)),1:length(folders)));
    if isreal(probeFolderIDs)
        
        for probe=1:length(probeFolderIDs)
            probeFolder=fullfile(sessionFolder,folders(probeFolderIDs(probe)).name);
            shankDir=dir(probeFolder);
            probeNeur=0;
            NP=NP+1;
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
            shankN=NaN(shanks,1);
            for shank=1:shanks
                
                shankFolder=fullfile(probeFolder,append('shank',num2str(shankList(shank))));
                clusterNames=readNPY(fullfile(shankFolder,'clusters.names.npy'));
                incClusters=clusterNames(readNPY(fullfile(shankFolder,'clusters.passedQualityControl.npy')));
                channelLocationsID = readNPY(fullfile(shankFolder,'channels.CCFregionIDs.npy'));
                channelXYZ = readNPY(fullfile(shankFolder,'channels.ML_AP_DV_fromBregma.npy'));
                channelLocations = textscan(fopen(fullfile(shankFolder,'channels.CCFregions.txt')),'%s');
                fclose('all');
                channelLocations = channelLocations{1};
                spikeTimes = readNPY(fullfile(shankFolder,'spikes.times.npy'));
                spikeClusters = readNPY(fullfile(shankFolder,'spikes.clusterNames.npy'));
                clusterChannel = readNPY(fullfile(shankFolder,'clusters.channels.npy'));
                probeNeur=probeNeur+length(incClusters);
                shankN(shank)=length(incClusters);
                
                %get clusters
                for clus=1:length(incClusters)
                    cluster = incClusters(clus);
                    regionName = channelLocations{clusterChannel(clusterNames==cluster)};
                    charOnly = regexp(regionName,'[c-zA-Z]','match');
                    noNums=cat(2,charOnly{:});
                    if sum(strcmp(allAreas,noNums(1:min([3 length(noNums)]))))==1
                        clusterTimes = double(spikeTimes(spikeClusters==cluster))/30000;
                        hz=length(clusterTimes)/(double(spikeTimes(end))/30000);
                        
                        %if hz>0.1
                        NN=NN+1;
                        neuronXYZ(NN,1:3)=channelXYZ(clusterChannel(clusterNames==cluster),:);
                        neuronRegion{NN,1}=channelLocations{clusterChannel(clusterNames==cluster)};
                        neuronRegionAdj{NN,1}=noNums(1:min([3 length(noNums)]));
                        neuronSubregion{NN,1}=noNums(1:min([5 length(noNums)]));
                        neuronGroup(NN,1)=0;
                        for grp=1:length(regionGroups)
                            if sum(strcmp(regionGroups{grp},noNums(1:min([3 length(noNums)]))))==1 neuronGroup(NN,1)=grp; end
                        end
                        
                        if strcmp(neuronRegionAdj{NN,1},'MOs')
                            almdist=sqrt((abs(neuronXYZ(NN,1)/1000)-1.5)^2+(neuronXYZ(NN,2)/1000-2.5)^2);
                            if almdist<=0.75 %1.5mm diameter, chen et al 2017
                                neuronRegionAdj{NN,1}='ALM';
                            end
                        end
                        
                        neuronHz(NN,1)=hz;
                        neuronSubject{NN,1}=sessionSubject{session};
                        neuronDate{NN,1}=sessionDate{session};
                        neuronCluster(NN,1)=cluster;
                        neuronSession(NN,1)=session;
                        neuronProbe(NN,1)=NP;
                        neuronRegionID(NN,1)=channelLocationsID(clusterChannel(clusterNames==cluster));
                        %end
                        
                        
                    end
                end
                
            end
        end
        probeNeurons(NP,1)=probeNeur;
    end
    sessionNeurons(session,1)=NN-NNstart;
    if NN>NNstart includedSessions(session,1)=true; end
    fprintf('Session %d \n',session);
end

totalNeurons=NN;

toc

%% analyze behavior
%run twice to see value lick plot

%parameters
binSize = 0.1;
binSizeSession=0.02;
window = [-2 8];
binEdges = window(1):binSize:window(2);
binTimes = window(1)+binSize/2:binSize:window(2)-binSize/2;

figure;

%parameters
ymax=6.5;
stimDuration = 1.5;
rewOnset = 2.5;
smoothing=8; %controls standard deviation of smoothing filter
plotrange = [-0.5 6];
%smoothing filter for licking PSTH
licksmoothbins=10; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',smoothing);
lickfilterweights=pdf(halfnormal,0:licksmoothbins);
subjects=unique(neuronSubject);

anticipatoryLicks={};
lickPSTH3={};
lickPSTH4={};
lickPSTHv={};
sessInd=find(includedSessions);
alldat=[];
allpreds=[];
allsubj=[];
trialsback=10;
for incses=1:sum(includedSessions)
    session=sessInd(incses);
    
    %get events
    sessionFolder=fullfile(direc,sessionSubject{session},sessionDate{session});
    lickTimes = readNPY(fullfile(sessionFolder,'licks.times.npy'));
    reward = readNPY(fullfile(sessionFolder,'rewards.times.npy'));
    cueTimes = readNPY(fullfile(sessionFolder,'cues.times.npy'));
    rewprob = readNPY(fullfile(sessionFolder,'cues.rewardProbability.npy'));
    odorset = readNPY(fullfile(sessionFolder,'cues.odorSet.npy')); 
    csp1 = cueTimes(rewprob==1&odorset==1);
    csp2 = cueTimes(rewprob==1&odorset==2);
    csf1 = cueTimes(rewprob==0.5&odorset==1);
    csf2 = cueTimes(rewprob==0.5&odorset==2);
    csm1 = cueTimes(rewprob==0&odorset==1);
    csm2 = cueTimes(rewprob==0&odorset==2);
    trialNo = [1:length(cueTimes)]';

    %get licks/bin
    binnedLicks=NaN(length(cueTimes),length(binEdges)-1);
    antLicks=NaN(length(cueTimes),1);
    conLicks=NaN(length(cueTimes),1);
    preLicks=NaN(length(cueTimes),1);
    for trial = 1:length(cueTimes)
        lickTimesRel = (lickTimes - cueTimes(trial));
        binnedLicks(trial,:)=histcounts(lickTimesRel,binEdges)/binSize;
        antLicks(trial,1) = sum(lickTimesRel>0 & lickTimesRel<2.5);
        conLicks(trial,1) = sum(lickTimesRel>3 & lickTimesRel<4.5);
        preLicks(trial,1) = sum(lickTimesRel>-2.5 & lickTimesRel<0);
    end
    
    trialTbl=table();
    trialTbl.trialNo=trialNo;
    trialTbl.csp1=ismember(round(cueTimes),round(csp1));
    trialTbl.csf1=ismember(round(cueTimes),round(csf1));
    trialTbl.csm1=ismember(round(cueTimes),round(csm1));
    trialTbl.csp2=ismember(round(cueTimes),round(csp2));
    trialTbl.csf2=ismember(round(cueTimes),round(csf2));
    trialTbl.csm2=ismember(round(cueTimes),round(csm2));    
    trialTbl.reward=ismember(round(cueTimes),round(reward-2.5));
    trialTbl.antLicks=antLicks;
    trialTbl.preLicks=preLicks;
    trialTbl.conLicks=conLicks;
    trialData{incses,1}=trialTbl;
    
    
    cuep1=trialTbl.csp1==1;
    cuep2=trialTbl.csp2==1;
    cuef1=trialTbl.csf1==1;
    cuef2=trialTbl.csf2==1;
    cue1=trialTbl.csp1==1 | trialTbl.csp2==1;
    cue2=trialTbl.csf1==1 | trialTbl.csf2==1;
    cue3=trialTbl.csm1==1 | trialTbl.csm2==1;
    
    dat=(trialTbl.antLicks)/max(trialTbl.antLicks);
    rew=trialTbl.reward;
    
    prevRew=zeros(length(rew),trialsback);
    prevCueRew=zeros(length(rew),trialsback);
    cuef1rew=rew(cuef1);
    cuef2rew=rew(cuef2);
    
    cuef1prev=ones(sum(cuef1),trialsback)/2;
    cuef2prev=ones(sum(cuef2),trialsback)/2;
    
    for tb=1:trialsback
        cuef1prev(tb+1:end,tb)=cuef1rew(1:end-tb);
        cuef2prev(tb+1:end,tb)=cuef2rew(1:end-tb);
        prevRew(tb+1:end,tb)=rew(1:end-tb);
    end
    prevCueRew(cue1,:)=1;
    prevCueRew(cuef1,:)=cuef1prev;
    prevCueRew(cuef2,:)=cuef2prev;
    p5r=sum(prevCueRew(:,1:5),2);
    %how many of each cue experienced
    num1=zeros(length(rew),1);
    num1(cue1)=[1:sum(cue1)]/sum(cue1);
    num2=zeros(length(rew),1);
    num2(cue2)=[1:sum(cue2)]/sum(cue2);
    num3=zeros(length(rew),1);
    num3(cue3)=[1:sum(cue3)]/sum(cue3);
    
    preds=[prevRew prevCueRew cue1 cue2 num1 num2 num3];
    
    
    subj=ones(length(rew),1)*find(ismember(subjects,sessionSubject(session)));
    allsubj=[allsubj;subj];
    alldat=[alldat;dat];
    allpreds=[allpreds;preds];

    %smooth lick traces
    smoothedLicks=[];
    for trial=1:length(cueTimes)
        for l=1:length(binTimes)
            smoothedLicks(trial,l)=sum(binnedLicks(trial,l-min([l-1 licksmoothbins]):l).*fliplr(lickfilterweights(1:min([l licksmoothbins+1]))))/sum(lickfilterweights(1:min([l licksmoothbins+1])));
        end
    end
    
    %3 conditions
    selections={};
    selections{1,1}=csp1; %cs+
    selections{2,1}=csf1; %cs50
    selections{3,1}=csm1; %cs-
    for condition=1:length(selections)
        sel = ismember(round(cueTimes),round(selections{condition}));
        lickPSTH3{condition,1}(incses,:)=nanmean(smoothedLicks(sel,:));
        anticipatoryLicks{condition,1}(incses,1)=nanmean(antLicks(sel)) - nanmean(preLicks(sel));    
    end
    selections{1,1}=csp2; %cs+
    selections{2,1}=csf2; %cs50
    selections{3,1}=csm2; %cs-
    for condition=1:length(selections)
        sel = ismember(round(cueTimes),round(selections{condition}));
        lickPSTH3{condition,2}(incses,:)=nanmean(smoothedLicks(sel,:));
        anticipatoryLicks{condition,1}(incses,2)=nanmean(antLicks(sel)) - nanmean(preLicks(sel));    
    end    
    
    %4 conditions
    selections{1,1}=csp1; %cs+
    selections{2,1}=csf1(ismember(round(csf1),round(reward-2.5))); %reward+
    selections{3,1}=csf1(~ismember(round(csf1),round(reward-2.5))); %reward-
    selections{4,1}=csm1; %cs-
    for condition=1:length(selections)
        sel = ismember(round(cueTimes),round(selections{condition}));
        lickPSTH4{condition,1}(incses,:)=nanmean(smoothedLicks(sel,:));
    end
    selections{1,1}=csp2; %cs+
    selections{2,1}=csf2(ismember(round(csf2),round(reward-2.5))); %reward+
    selections{3,1}=csf2(~ismember(round(csf2),round(reward-2.5))); %reward-
    selections{4,1}=csm2; %cs-
    for condition=1:length(selections)
        sel = ismember(round(cueTimes),round(selections{condition}));
        lickPSTH4{condition,2}(incses,:)=nanmean(smoothedLicks(sel,:));
    end    
    
    %lick PSTH based on trial value from coefficient weights
    if exist('trialValues')
        %previous rewards
        cue1 = ismember(round(cueTimes),round([csp1;csp2]));
        cue2 = ismember(round(cueTimes),round([csf1;csf2]));
        selections{1,1}=trialValues{incses}<0.05; %CS- low
        selections{2,1}=trialValues{incses}>=0.05 & trialValues{incses}<0.1; %CS- high
        selections{3,1}=trialValues{incses}>=0.2 & trialValues{incses}<0.25; %CS50 low
        selections{4,1}=trialValues{incses}>=0.25 & trialValues{incses}<0.3;
        selections{5,1}=trialValues{incses}>=0.3 & trialValues{incses}<0.35;
        selections{6,1}=trialValues{incses}>=0.35 & trialValues{incses}<0.4;
        selections{7,1}=trialValues{incses}>=0.4 & trialValues{incses}<0.45;
        selections{8,1}=trialValues{incses}>=0.45 & cue2; %CS50 high
        selections{9,1}=trialValues{incses}<0.5 & cue1; %CS+ low
        selections{10,1}=trialValues{incses}>=0.5; %CS+ high
        for condition=1:length(selections)
            sel = selections{condition};
            lickPSTHv{condition,1}(incses,:)=nanmean(smoothedLicks(sel,:));
        end
    end
    
    
end

for os=1:2
%plot PSTHs
subplot(3,3,os)
hold on

colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];

xvals = find(binTimes>=plotrange(1) & binTimes<=2.5);
for condition=1:3
    
    psthsubj=NaN(length(subjects),size(lickPSTH3{condition,os},2));
    for subject=1:length(subjects)
        psthsubj(subject,:)=nanmean(lickPSTH3{condition,os}(ismember(sessionSubject,subjects(subject)),:));
    end
    
    %get values
    psth=nanmean(psthsubj);
    sem=nanste(psthsubj,1); %calculate standard error of the mean
    up=psth+sem;
    down=psth-sem;
    
    %plotting
    plot(binTimes(xvals),psth(xvals),'Color',colors{condition,1},'linewidth',1);
    patch([binTimes(xvals),binTimes(xvals(end):-1:xvals(1))],[up(xvals),down(xvals(end):-1:xvals(1))],colors{condition,1},'EdgeColor','none');alpha(0.5);
end


colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.6 0.1 0.6];
colors{3,1}=[0.25 0 0.25];
colors{4,1}=[0.3 0.3 0.3];

xvals = find(binTimes>=2.5 & binTimes<=plotrange(2));
p={};
for condition=1:4
    
    psthsubj=NaN(length(subjects),size(lickPSTH4{condition,os},2));
    for subject=1:length(subjects)
        psthsubj(subject,:)=nanmean(lickPSTH4{condition,os}(ismember(sessionSubject,subjects(subject)),:));
    end
    
    %get values
    psth=nanmean(psthsubj);
    sem=nanste(psthsubj,1); %calculate standard error of the mean
    up=psth+sem;
    down=psth-sem;
    
    %plotting
    p{condition}=plot(binTimes(xvals),psth(xvals),'Color',colors{condition,1},'linewidth',1);
    patch([binTimes(xvals),binTimes(xvals(end):-1:xvals(1))],[up(xvals),down(xvals(end):-1:xvals(1))],colors{condition,1},'EdgeColor','none');alpha(0.3);
end

patch([0 0 stimDuration stimDuration],[0 ymax ymax 0],[0.6 0.3 0],'edgecolor','none');alpha(0.3);
plot([rewOnset rewOnset],[0 ymax],'color',[0 0.6 0.3]);
axis([plotrange 0 ymax]);
xlabel('seconds from odor onset');
ylabel('Licks/s');

if os==2 legend([p{:}],'CS+','CS50(r)','CS50(u)','CS-','location','NW'); end
title(['Odor set ' num2str(os)]);
end

if exist('trialValues')
    colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];

    subplot(3,3,7);
    hold on;
    xvals = find(binTimes>=plotrange(1) & binTimes<=2.5);
    p={};
    for condition=1:3
        
        psthsubj=NaN(length(subjects),size(lickPSTH3{condition,1},2));
        for subject=1:length(subjects)
            psthsubj(subject,:)=nanmean(nanmean(cat(3,lickPSTH3{condition,1}(ismember(sessionSubject,subjects(subject)),:),lickPSTH3{condition,2}(ismember(sessionSubject,subjects(subject)),:)),3));
        end
        
        %get values
        psth=nanmean(psthsubj);
        sem=nanste(psthsubj,1); %calculate standard error of the mean
        up=psth+sem;
        down=psth-sem;
        
        %plotting
        p{condition}=plot(binTimes(xvals),psth(xvals),'Color',colors{condition,1},'linewidth',1.25);        
        patch([binTimes(xvals),binTimes(xvals(end):-1:xvals(1))],[up(xvals),down(xvals(end):-1:xvals(1))],colors{condition,1},'EdgeColor','none');alpha(0.05);
    end
    
    
    %patch([0 0 stimDuration stimDuration],[0 ymax ymax 0],[0.6 0.3 0],'edgecolor','none');alpha(0.3);
    axis([0 2.5 0 4.5]);
    xlabel('seconds from odor onset');
    ylabel('Licks/s');
    
    %legend([p{:}],'0','1','2','3','4','5');
    title('Cue value');
    
    
    %with trial history
    subplot(3,3,8);
    hold on;
valcolors{1,1}=[0.2 0.2 0.2];
valcolors{2,1}=[0.5 0.5 0.5];
valcolors{3,1}=[0.35 0.05 0.35];
valcolors{4,1}=[0.4 0.1 0.4];
valcolors{5,1}=[0.5 0.15 0.5];
valcolors{6,1}=[0.6 0.2 0.6];
valcolors{7,1}=[0.7 0.25 0.7];
valcolors{8,1}=[0.8 0.3 0.8];
valcolors{9,1}=[0.1 0.6 0.2];
valcolors{10,1}=[0.2 0.8 0.3];
    xvals = find(binTimes>=plotrange(1) & binTimes<=2.5);
    p={};
    for condition=1:10
        
        psthsubj=NaN(length(subjects),size(lickPSTHv{condition,1},2));
        for subject=1:length(subjects)
            psthsubj(subject,:)=nanmean(lickPSTHv{condition,1}(ismember(sessionSubject,subjects(subject)),:));
        end
        
        %get values
        psth=nanmean(psthsubj);
        sem=nanste(psthsubj,1); %calculate standard error of the mean
        up=psth+sem;
        down=psth-sem;
        
        %plotting
        patch([binTimes(xvals),binTimes(xvals(end):-1:xvals(1))],[up(xvals),down(xvals(end):-1:xvals(1))],valcolors{condition,1},'EdgeColor','none');alpha(0.05);
    end

    for condition=1:10
        
        psthsubj=NaN(length(subjects),size(lickPSTHv{condition,1},2));
        for subject=1:length(subjects)
            psthsubj(subject,:)=nanmean(lickPSTHv{condition,1}(ismember(sessionSubject,subjects(subject)),:));
        end
        
        %get values
        psth=nanmean(psthsubj);
        sem=nanste(psthsubj,1); %calculate standard error of the mean
        up=psth+sem;
        down=psth-sem;
        
        %plotting
        p{condition}=plot(binTimes(xvals),psth(xvals),'Color',valcolors{condition,1},'linewidth',1.25);
    end    
    
    
    %patch([0 0 stimDuration stimDuration],[0 ymax ymax 0],[0.6 0.3 0],'edgecolor','none');alpha(0.3);
    axis([0 2.5 0 4.5]);
    xlabel('seconds from odor onset');
    ylabel('Licks/s');
    
    %legend([p{:}],'0','1','2','3','4','5');
    title('Cue value');
end


colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];

subplot(3,6,6);
hold on;
mat4anov=[];
for cue=1:3
    antLicksSubj=NaN(length(subjects),2);
    for os=1:2
        for subject=1:length(subjects)
            antLicksSubj(subject,os)=nanmean(anticipatoryLicks{cue}(ismember(sessionSubject,subjects(subject)),os));
        end
    end
    errorbar([1 2],mean(antLicksSubj),nanste(antLicksSubj,1),'linewidth',1.5,'color',colors{cue});
    plot([1 2],antLicksSubj,'color',colors{cue},'linewidth',0.5);
    mat4anov=cat(1,mat4anov,antLicksSubj);
end
ylabel('\Delta anticipatory licks');
plot([0 3],[0 0],':','color','k','linewidth',0.75);
xticks([1 2]);
xticklabels({'1','2'});
xlabel('odor set');
xlim([0.75 2.25]);
%[p,tbl,stats]=anova2(mat4anov,length(subjects));
%c=multcompare(stats,'estimate','row');

rweights=[];
cweights=[];
for s=1:length(subjects)
sel=allsubj==s;
mdl=fitlm(allpreds(sel,:),alldat(sel));
rweights(s,:)=mdl.Coefficients.Estimate(2:trialsback+1);
cweights(s,:)=mdl.Coefficients.Estimate(trialsback+2:trialsback*2+1);
end

%allmice together
sel=allsubj>0;
mdl=fitlm(allpreds(sel,:),alldat(sel));
subplot(3,3,4)
hold on;
errorbar(1:trialsback,mdl.Coefficients.Estimate(2:trialsback+1),mdl.Coefficients.SE(2:trialsback+1),'linewidth',1.5,'color','k');
plot(1:trialsback,rweights,'linewidth',0.5);
plot([0 11],[0 0],':','color','k');
ylim([-0.05 0.15]);
xlim([0 trialsback+1]);
xticks([1 5 10]);
ylabel('coefficient');
xlabel('trials back');
title('previous trial rewarded');

subplot(3,3,5)
hold on;
errorbar(1:trialsback,mdl.Coefficients.Estimate(trialsback+2:trialsback*2+1),mdl.Coefficients.SE(trialsback+2:trialsback*2+1),'linewidth',1.5,'color','k');
plot(1:trialsback,cweights,'linewidth',0.5);
plot([0 11],[0 0],':','color','k');
ylim([-0.05 0.15]);
xlim([0 trialsback+1]);
xticks([1 5 10]);

ylabel('coefficient');
xlabel('trials back');
title('previous cue rewarded');

if ~exist('trialValues')
%get values from weights
trialValues={};
for session=1:sum(includedSessions)
    
    
    trialTbl=trialData{session};
    cuep1=trialTbl.csp1==1;
    cuep2=trialTbl.csp2==1;
    cuef1=trialTbl.csf1==1;
    cuef2=trialTbl.csf2==1;
    cue1=trialTbl.csp1==1 | trialTbl.csp2==1;
    cue2=trialTbl.csf1==1 | trialTbl.csf2==1;
    cue3=trialTbl.csm1==1 | trialTbl.csm2==1;
    
    dat=(trialTbl.antLicks)/max(trialTbl.antLicks);
    rew=trialTbl.reward;
    
    prevRew=zeros(length(rew),trialsback);
    prevCueRew=zeros(length(rew),trialsback);
    cuef1rew=rew(cuef1);
    cuef2rew=rew(cuef2);
    
    cuef1prev=ones(sum(cuef1),trialsback)/2;
    cuef2prev=ones(sum(cuef2),trialsback)/2;
  
    for tb=1:trialsback
        cuef1prev(tb+1:end,tb)=cuef1rew(1:end-tb);
        cuef2prev(tb+1:end,tb)=cuef2rew(1:end-tb);
        prevRew(tb+1:end,tb)=rew(1:end-tb);
    end
    prevCueRew(cue1,:)=1;
    prevCueRew(cuef1,:)=cuef1prev;
    prevCueRew(cuef2,:)=cuef2prev;
    
    preds=[prevRew prevCueRew cue1 cue2];
    trialValues{session,1}=preds * mdl.Coefficients.Estimate(2:end-3);
    

end
end

%% do analysis on all neurons
preRun=true;
analysisFile='value20220610.mat';

analysisWindow=[-1 2.5];
binsPerTrial=diff(analysisWindow)/binSize;

%for shifting lick vector in time
%each shift is binSize*2
numShiftsB=3; %making licks happen later than reality
numShiftsF=2; %making licks happen earlier than reality

%base model
predictornames{1}='csp1';
predictornames{2}='csf1';
predictornames{3}='csm1';
predictornames{4}='csp2';
predictornames{5}='csf2';
predictornames{6}='csm2';
predictornames{7}='licks';
predictornames{8}='boutStart';
predictornames{9}='lickRate';
%windows for these predictors
windows{1}=[0 2.5];
windows{2}=[0 2.5];
windows{3}=[0 2.5];
windows{4}=[0 2.5];
windows{5}=[0 2.5];
windows{6}=[0 2.5];
windows{7}=[-0.3 0.3];
windows{8}=[-0.3 2];
windows{9}=[-numShiftsB*binSize (numShiftsF+1)*binSize]; %there are numShifts*2+1 total for this, so fudging it a bit
windows{10}=[0 binSize*6]; %block constants

%model with just value for cue
predictornames{1}='cueValue';
predictornames{2}='licks';
predictornames{3}='boutStart';
predictornames{4}='lickRate';
%windows for these predictors
windows3{1}=[0 2.5];
windows3{2}=[-0.3 0.3];
windows3{3}=[-0.3 2];
windows3{4}=[-numShiftsB*binSize (numShiftsF+1)*binSize]; %there are numShifts*2+1 total for this, so fudging it a bit
windows3{5}=[0 binSize*6]; %block constants

binPred=[];
binsofar=0;
for win=1:length(windows)
    winbins=diff(windows{win})/binSize;
    binPred(binsofar+1:binsofar+winbins)=win;
    binsofar=binsofar+winbins;
end
lickPreds=ismember(binPred,[7 8 9]);
cuePreds=ismember(binPred,[1 2 3 4 5 6]);
conPreds=ismember(binPred,10);

submodelsels={};
submodelsels{1}=cuePreds==0;
submodelsels{2}=lickPreds==0;

tic

if preRun
    load(analysisFile)
else

Xall=cell(sum(includedSessions),1);
Xall3=cell(sum(includedSessions),1);
Yall=cell(sum(includedSessions),1);
YallSh=cell(sum(includedSessions),1);
PSTH3=cell(3,2);
PSTH4=cell(4,2);
PSTH6=cell(6,2);
PSTHv=cell(10,1);

NS=0;

%smoothing filter for licks
smoothbinsLick=25; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',20);
filterweightsLick=pdf(halfnormal,0:smoothbinsLick);

%smoothing filter for individual trials
smoothbinsTrl=50; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',15);
filterweightsTrl=pdf(halfnormal,0:smoothbinsTrl);


NN=0;
sessInd=find(includedSessions);

for incses=1:sum(includedSessions)
    session=sessInd(incses);
    
    NS=NS+1;
    %get events
    sessionFolder=fullfile(direc,sessionSubject{session},sessionDate{session});
    lickTimes = readNPY(fullfile(sessionFolder,'licks.times.npy'));
    reward = readNPY(fullfile(sessionFolder,'rewards.times.npy'));
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
    
    %first lick in bout
    ibi=0.5;
    boutStart = lickTimes(cat(1,1,1+find(diff(lickTimes)>ibi)));
    
    %first rewarded lick
    rewardswithlicks=reward(reward<max(lickTimes));
    lickreward = arrayfun(@(x) lickTimes(find(lickTimes>reward(x),1)),1:length(rewardswithlicks))';
    lickreward = unique(lickreward);
    
    %only use lick reward when at least 400ms after reward
    lrdelay=[];
    for l=1:length(lickreward)
        ldiffs=lickreward(l)-reward;
        lrdelay(l,1)=min(ldiffs(ldiffs>0));
    end  
    lickreward(lrdelay<0.4)=[];
    
    rewarded = arrayfun(@(x) cueTimes(find(cueTimes<reward(x),1,'last')),1:length(reward))';
    rewarded = ismember(cueTimes,rewarded);


    
    binEdgesSession=cueTimes(1)-5:binSizeSession:cueTimes(end)+10;
    binTimesSession=cueTimes(1)-5+binSizeSession/2:binSizeSession:cueTimes(end)+10-binSizeSession/2;
    binnedLicksGLM=histcounts(lickTimes,binEdgesSession) / binSizeSession;
    lickRateGLM=NaN(length(binnedLicksGLM),1);
    for l=1:length(binnedLicksGLM)
        lickRateGLM(l,1)=sum(binnedLicksGLM(l-min([l-1 smoothbinsLick]):l).*fliplr(filterweightsLick(1:min([l smoothbinsLick+1]))))/sum(filterweightsLick(1:min([l smoothbinsLick+1])));
    end

    binnedLicksGLM=NaN(length(cueTimes),length(binEdges)-1);
    binTimesGLM=NaN(length(cueTimes),length(binEdges)-1);
    for trial = 1:length(cueTimes)
        binTimesRel = (binTimesSession - cueTimes(trial));
        for bin=1:length(binTimes)
            binnedLicksGLM(trial,bin)=mean(lickRateGLM(binTimesRel>=binEdges(bin)&binTimesRel<binEdges(bin+1)));
            binTimesGLM(trial,bin)=binEdges(bin)+cueTimes(trial)+binSize/2;
        end
    end
    binnedLicksGLM(isnan(binnedLicksGLM))=0;  %this only happens if bins go beyond end of session, so extremely rare

    vector1={};
    includedBins=find(binTimes>analysisWindow(1) & binTimes<analysisWindow(2));
    for shift=1:numShiftsB+numShiftsF+1

        includedBinsShifted=includedBins+shift*2-numShiftsB*2-2;
        includedLickBins=(binnedLicksGLM(:,includedBinsShifted))';
        vector1{shift}=includedLickBins(:)/max(includedLickBins(:));
    end
    vector1=vector1(end:-1:1); %flipping lick vector order so first one is the "earliest"
    binTimesGLMincluded=binTimesGLM(:,includedBins)';
    GLMbins{NS,1} = binTimesGLMincluded(:);

    %constant for each block
    ITIs=diff(cueTimes);
    [ITIsort,order]=sort(ITIs);
    biggestITIs=order(end-4:end);
    lastTrialBlock=sort(biggestITIs);
    lastTrialBlock=[0;lastTrialBlock;length(cueTimes)];
    blockSize=diff(lastTrialBlock);
    
    blockCs={};
    for block=1:6
        blockC=zeros(length(vector1{1}),1);
        blockC(binsPerTrial*lastTrialBlock(block)+1:binsPerTrial*lastTrialBlock(block)+binsPerTrial*blockSize)=1;
        blockCs{block}=blockC;
    end
    
    %make predictor matrix (called A here)
    continuousPredictors=[vector1 blockCs];
    discreteUnfixed.times={csp1,csf1,csm1,csp2,csf2,csm2,lickTimes,boutStart};
    discreteUnfixed.windows=windows(1:8);
    discreteUnfixed.binTimes=binTimesGLMincluded(:);
    discreteUnfixed.values={ones(length(csp1),1),ones(length(csf1),1),ones(length(csm1),1),...
        ones(length(csp2),1),ones(length(csf2),1),ones(length(csm2),1),...
        ones(length(lickTimes),1),ones(length(boutStart),1)};
    A=makeA(discreteUnfixed,continuousPredictors);
    Xall{NS,1}=A;
    
    %make predictor matrix for value cue
    discreteUnfixed.times={cueTimes,lickTimes,boutStart};
    discreteUnfixed.windows=windows3(1:3);
    discreteUnfixed.binTimes=binTimesGLMincluded(:);
    discreteUnfixed.values={trialValues{incses},ones(length(lickTimes),1),ones(length(boutStart),1)};
    A3=makeA(discreteUnfixed,continuousPredictors);
    Xall3{NS,1}=A3;
    
    folders=dir(sessionFolder);
    probeFolderIDs = find(arrayfun(@(x) strcmp('p',folders(x).name(1)),1:length(folders)));
    if isreal(probeFolderIDs)
        NNses=0;
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
                        NNses=NNses+1;
                        NN=NN+1;
                        clusterTimes = double(spikeTimes(spikeClusters==cluster))/30000;
                        
                        
                        %get spike trace
                        spikeRate = histcounts(clusterTimes,binEdgesSession)/binSizeSession;
                        
                        %smooth
                        spikeRateSmooth=NaN(1,length(spikeRate));
                        for l=1:length(spikeRateSmooth)
                            spikeRateSmooth(1,l)=sum(spikeRate(1,l-min([l-1 smoothbinsTrl]):l).*fliplr(filterweightsTrl(1:min([l smoothbinsTrl+1]))))/sum(filterweightsTrl(1:min([l smoothbinsTrl+1])));
                        end
                        
                        
                        binnedActivity=NaN(length(cueTimes),length(binEdges)-1);
                        numBins=length(binTimes);
                        for trial = 1:length(cueTimes)
                            binTimesRel = (binTimesSession - cueTimes(trial));
                            [binsPerBin,~,binNumber]=histcounts(binTimesRel,binEdges);
                            binsUsed=binNumber>0;
                            binnedActivity(trial,:)=accumarray(binNumber(binsUsed)',spikeRateSmooth(binsUsed)',[numBins 1]) ./ binsPerBin';
                        end
                        binnedActivity(isnan(binnedActivity))=0;binnedActivity(binnedActivity==inf)=0; %this only happens if bins go beyond end of session, so extremely rare
                        
                        %to avoid silly z-scores
                        neuronStd(NN,1)=std(binnedActivity(:,binTimes<0),0,'all');
                        if std(binnedActivity(:,binTimes<0),0,'all')>=1
                            normActivity=(binnedActivity-mean(binnedActivity(:,binTimes<0),'all')) / std(binnedActivity(:,binTimes<0),0,'all');
                        else
                            normActivity=(binnedActivity-mean(binnedActivity(:,binTimes<0),'all'));
                        end
                        %3 conditions
                        selections={};
                        selections{1,1}=csp1; %cs+
                        selections{2,1}=csf1; %cs50
                        selections{3,1}=csm1; %cs-
                        for condition=1:length(selections)
                            sel = ismember(round(cueTimes),round(selections{condition}));
                            PSTH3{condition,1}(NN,:)=nanmean(normActivity(sel,:));
                        end
                        selections{1,1}=csp2; %cs+
                        selections{2,1}=csf2; %cs50
                        selections{3,1}=csm2; %cs-
                        for condition=1:length(selections)
                            sel = ismember(round(cueTimes),round(selections{condition}));
                            PSTH3{condition,2}(NN,:)=nanmean(normActivity(sel,:));
                        end
                        
                        %4 conditions
                        selections={};
                        selections{1,1}=csp1; %cs+
                        selections{2,1}=csf1(ismember(round(csf1),round(reward-2.5))); %reward+
                        selections{3,1}=csf1(~ismember(round(csf1),round(reward-2.5))); %reward-
                        selections{4,1}=csm1; %cs-
                        for condition=1:length(selections)
                            sel = ismember(round(cueTimes),round(selections{condition}));
                            PSTH4{condition,1}(NN,:)=nanmean(normActivity(sel,:));
                        end
                        selections{1,1}=csp2; %cs+
                        selections{2,1}=csf2(ismember(round(csf2),round(reward-2.5))); %reward+
                        selections{3,1}=csf2(~ismember(round(csf2),round(reward-2.5))); %reward-
                        selections{4,1}=csm2; %cs-
                        for condition=1:length(selections)
                            sel = ismember(round(cueTimes),round(selections{condition}));
                            PSTH4{condition,2}(NN,:)=nanmean(normActivity(sel,:));
                        end
                        
                        %3 conditions, odd and even trials
                        selections={};
                        selections{1,1}=csp1; %cs+
                        selections{2,1}=csp1; %cs+
                        selections{3,1}=csf1; %cs50
                        selections{4,1}=csf1; %cs50
                        selections{5,1}=csm1; %cs-
                        selections{6,1}=csm1; %cs-
                        for condition=1:length(selections)
                            sel = ismember(round(cueTimes),round(selections{condition})) & rem(trialNo',2)==rem(condition,2);
                            PSTH6{condition,1}(NN,:)=nanmean(normActivity(sel,:));
                        end
                        selections{1,1}=csp2; %cs+
                        selections{2,1}=csp2; %cs+
                        selections{3,1}=csf2; %cs50
                        selections{4,1}=csf2; %cs50
                        selections{5,1}=csm2; %cs-
                        selections{6,1}=csm2; %cs-
                        for condition=1:length(selections)
                            sel = ismember(round(cueTimes),round(selections{condition})) & rem(trialNo',2)==rem(condition,2);
                            PSTH6{condition,2}(NN,:)=nanmean(normActivity(sel,:));
                        end
                        
                        %value
                        cue1 = ismember(round(cueTimes),round([csp1;csp2]));
                        cue2 = ismember(round(cueTimes),round([csf1;csf2]));
                        selections{1,1}=trialValues{incses}<0.05; %CS- low
                        selections{2,1}=trialValues{incses}>=0.05 & trialValues{incses}<0.1; %CS- high
                        selections{3,1}=trialValues{incses}>=0.2 & trialValues{incses}<0.25; %CS50 low
                        selections{4,1}=trialValues{incses}>=0.25 & trialValues{incses}<0.3;
                        selections{5,1}=trialValues{incses}>=0.3 & trialValues{incses}<0.35;
                        selections{6,1}=trialValues{incses}>=0.35 & trialValues{incses}<0.4;
                        selections{7,1}=trialValues{incses}>=0.4 & trialValues{incses}<0.45;
                        selections{8,1}=trialValues{incses}>=0.45 & cue2; %CS50 high
                        selections{9,1}=trialValues{incses}<0.5 & cue1; %CS+ low
                        selections{10,1}=trialValues{incses}>=0.5; %CS+ high
                        for condition=1:length(selections)
                            sel = selections{condition};
                            PSTHv{condition,1}(NN,:)=nanmean(normActivity(sel,:));
                        end
                        
                        
                        %%% prepare for GLM %%%%
                        includedActivity=normActivity(:,includedBins)';
                        shOrder=randperm(size(normActivity,1));
                        includedActivitySh=normActivity(shOrder,includedBins)';
                        activityVector=includedActivity(:);
                        activityVectorSh=includedActivitySh(:);
                        
                        Yall{NS,1}(:,NNses)=activityVector;
                        YallSh{NS,1}(:,NNses)=activityVectorSh;
                        
                    end
                end
                
                
            end
        end
    end
    
    fprintf('Session %d \n',session);
end

save('newAnalysis.mat','PSTH3','PSTH4','PSTH6','PSTHv','PSTHf','Xall','Xall2','Xall3','Yall','YallSh','-v7.3');

end

toc

%% get reduced rank predictors

[~, bR, R2] = CanonCor2all(Yall, Xall);
[~, bR3, R2] = CanonCor2all(Yall, Xall3);

%% perform regression with reduced rank predictors
preRun=true;
glmFile='valueGLM20230127';

tic

if preRun
    load(glmFile);
else

components=10;
folds=4;
NS=0;
varExp=NaN(totalNeurons,7);
varExpValueSh=NaN(totalNeurons,94);
predF=cell(totalNeurons,1);
kernels=NaN(totalNeurons,1);
lambdaVal=NaN(totalNeurons,1);
lambda=[0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.5];
%lambda=[0.001 0.01 0.1 0.5];
binsPerCue=diff(windows{1})/binSize;
NN=0;

for session=1:sum(includedSessions)

        NS=NS+1;
        
        A = Xall{NS};
        rA = A * bR(:,1:components);
        trains = getfolds(rA,folds,binsPerTrial);
        
        %submodels
        rAs={};
        As={};
        for sub=1:length(submodelsels)
            At = A;
            At(:,submodelsels{sub}==0)=0;
            As{sub} = At;
            rAs{sub} = At * bR(:,1:components);
        end
            
        varExpLam=[];
        for neuron=1:size(Yall{NS},2)
            NN=NN+1;
            y=Yall{NS}(:,neuron);
            
            %cross-validated variance explained
            predLam=NaN(size(y,1),length(lambda));
            for fold=1:folds
                train=trains{fold};
                test=train==0;         
                fitK=lassoglm(rA(train,:),y(train),'normal','alpha',0.5,'lambda',lambda);
                kernel=bR(:,1:components)*fitK;
                predLam(test,:)=A(test,:)*kernel;
            end
            
            %get best lambda
            for lam=1:length(lambda)
            varExpLam(1,lam)=1- var(y-predLam(:,lam))/var(y);
            end
            [varExp(NN,1),ind] = max(varExpLam);
            thisLambda=lambda(ind);
            lambdaVal(NN,1)=thisLambda;
            
            
            pred=predLam(:,ind);
            if rem(NN,10)==0 predF{NN,1}=pred; end
            
  
            %full data to get kernels
            fitK=lassoglm(rA,y,'normal','alpha',0.5,'lambda',thisLambda);
            kernel=bR(:,1:components)*fitK;
            kernels(NN,1:length(kernel))=kernel;
            
            %submodels to find unique variance and total variance for each
            %variable
            for sub=1:length(submodelsels)
                pred=NaN(size(y,1),1);
                for fold=1:folds
                    train=trains{fold};
                    test=train==0;                  
                    fitK=lassoglm(rAs{sub}(train,:),y(train),'normal','alpha',0.5,'lambda',thisLambda);
                    kernel=bR(:,1:components)*fitK;
                    pred(test)=As{sub}(test,:)*kernel;
                end
                if rem(NN,10)==0 predF{NN,1+sub}=pred; end
                varExp(NN,1+sub) = 1- var(y-pred)/var(y);
                
            end
                      
            
        end
        fprintf('Session #%d \n',NS);
    
end

save('newGLMFile.mat','varExp','kernels','lambdaVal','predF','-v7.3');

end
toc


%% perform value analysis with a ton of shuffles
runglm=false;
%generate all value shuffles
valueShuffles = [];
valueVersion = [];
valueLevels = [1 1 0.7 0.7 0.1 0.1;1 1 0.5 0.5 0 0;1 1 1 1 0 0;1 1 1 0 0 0;1 1 0 0 0 0;1 0 0 0 0 0;1 1 1 1 1 0;1 1 1 1 1 1];
for vl=1:size(valueLevels,1)
allCombos = perms(valueLevels(vl,:));
allCombos = unique(allCombos,'rows');
valueShuffles = cat(1,valueShuffles,allCombos);
valueVersion = cat(1,valueVersion,vl*ones(size(allCombos,1),1));
end
valueShuffles=fliplr(valueShuffles);

folds=4;
NNstart=0;
%varExpValueNew=NaN(totalNeurons,length(valueShuffles));
%kernels=NaN(totalNeurons,31);
binsPerCue=diff(windows{1})/binSize;

opts.alpha=0.5;
glmopts=glmnetSet(opts);

if runglm

sessInd=find(includedSessions);
tic
for NS=1:sum(includedSessions)
    session=sessInd(incses);
    sessionFolder=fullfile(direc,sessionSubject{session},sessionDate{session});
    cueTimes = readNPY(fullfile(sessionFolder,'cues.times.npy'));
    rewprob = readNPY(fullfile(sessionFolder,'cues.rewardProbability.npy'));
    odorset = readNPY(fullfile(sessionFolder,'cues.odorSet.npy'));
    cueOrder=zeros(length(cueTimes),1);
    cueOrder(rewprob==1&odorset==1)=1;
    cueOrder(rewprob==1&odorset==2)=2;
    cueOrder(rewprob==0.5&odorset==1)=3;
    cueOrder(rewprob==0.5&odorset==2)=4;
    cueOrder(rewprob==0&odorset==1)=5;
    cueOrder(rewprob==0&odorset==2)=6;
    
    A = Xall{NS}(:,[1:binsPerCue size(Xall{NS},2)-5:size(Xall{NS},2)]);
    trains = getfolds(A,folds,binsPerTrial);
    for s=length(valueShuffles)-6:length(valueShuffles) %1:length(valueShuffles)
        for c=1:length(cueTimes)
            cueValue=valueShuffles(s,cueOrder(c));
            dm=cueValue*eye(binsPerCue);
            A(c*binsPerTrial-binsPerCue+1:c*binsPerTrial,1:binsPerCue)=dm;
        end
        NN=NNstart;
        for neuron=1:size(Yall{NS},2)
            NN=NN+1;
            y=Yall{NS}(:,neuron);
            thisLambda=lambdaVal(NN,1);
            %cross-validated variance explained
            pred=NaN(size(y,1),1);
            for fold=1:folds
                train=trains{fold};
                test=train==0;
                %kernel=lassoglm(A(train,:),y(train),'normal','alpha',0.5,'lambda',thisLambda);
                fit = glmnet(A(train,:),y(train),[],glmopts);
                prediction = glmnetPredict(fit,A(test,:));
                pred(test)=prediction(:,end);
            end
            varExpValueNew(NN,s) = 1- var(y-pred)/var(y);
            kernel = fit.beta(:,end);
            if s==1 kernels(NN,:)=kernel; end
        end           
    end
    NNstart=NN;
    fprintf('Session #%d \n',NS);
end
toc
save('valueGLMFile.mat','varExp','varExpValueNew','kernels','lambdaVal','predF','-v7.3');

end

%% value coding analysis
improvCutoff=0.02;
overallCutoff=0.02;
improvement=varExp(:,1)-varExp;

%get proportions
predictive = varExp(:,1)>overallCutoff; %more than 2% of variance explained
behavior = improvement(:,3)>improvCutoff; %num licks or lick rate
cueK = improvement(:,2)>improvCutoff; %if best model is at least cutoff better than no cue kernels

totalCueVariance=varExp(:,1)-varExp(:,2); %best model compared to no cue model
valueImp=varExp(:,1)-varExpValueNew; %how much better the best model is than value, shuffles, and one kernel
%correct it so you can't be worse than no cues at all (should I do this?)
for n=1:length(valueImp)
    valueImp(n,valueImp(n,:)>totalCueVariance(n))=totalCueVariance(n);
end

includedModels = [91:length(valueShuffles)];
%includedModels = [1 91 2:90 92:length(valueShuffles)];
%includedModels = [1:90 181:length(valueShuffles)];
%includedModels = [1:length(valueShuffles)];
shuffleLabels = cell(length(includedModels),1);
for sl=1:length(includedModels)
    m=includedModels(sl);
    shuffleLabels{sl} = sprintf('%g,%g,%g,%g,%g,%g',valueShuffles(m,1),valueShuffles(m,2),valueShuffles(m,3),valueShuffles(m,4),valueShuffles(m,5),valueShuffles(m,6));
end
[~,bestModel]=max(varExpValueNew(:,includedModels),[],2);
% [~,bestModel]=max(varExpValueNew,[],2);
% bestModel(bestModel>90)=bestModel(bestModel>90)-90;
totalModels = length(includedModels);

figure;

corrWVal=corrcoef(valueShuffles(includedModels,:)');
corrWVal(isnan(corrWVal))=0;
[valCorr,valRank] = sort(corrWVal(:,1),'descend');
[~,mdlPosition] = sort(valRank);

subplot(5,1,1);
hold on
frequency=histcounts(bestModel(cueK&~behavior&~isnan(bestModel)),0.5:1:totalModels+0.5,'normalization','probability');
%b=bar(1:max(bestModel),frequency,'facecolor',[0 0.2 0.4],'linewidth',0.01);
b=bar(mdlPosition,frequency,'facecolor',[0 0.2 0.4],'linewidth',0.01);

b.FaceColor='flat';
for e=2:17
%b.CData(valRank(e),:)=[0.7 0 1]; %trial type
b.CData(e,:)=[0.7 0 1]; %trial type
end
b.CData(1,:)=[0 0.7 1]; %value
b.CData(mdlPosition(end),:)=[0.6 0.6 0.6]; %non-specific
%b.CData(end,:)=[0.6 0.6 0.6]; %non-specific
ylabel('fraction of cue cells');

chance=1/totalModels;
plot([0 totalModels+1],[chance chance],':','color','k','linewidth',0.5);
xlim([0 totalModels+1]);
ylim([0 0.2]);
%xticks(mdlPosition(topModels(1:end-1)));
xticks([]);
yticks(0:0.1:0.3);
xlabel('cue coding models, sorted by similarity to value model');

subplot(4,5,6);
hold on
scatter(valCorr(2:17,1),frequency(valRank(2:17)),24,[0.7 0 1],'filled');
scatter(valCorr(1,1),frequency(valRank(1)),24,[0 0.7 1],'filled');
scatter(valCorr(18:end,1),frequency(valRank(18:end)),24,[0 0.2 0.4],'filled');
scatter(valCorr(mdlPosition(end),1),frequency(end),24,[0.6 0.6 0.6],'filled');
plot([-1 1],[1/length(includedModels) 1/length(includedModels)],':','color','k','linewidth',0.5);
ylim([0 0.2]);
ylabel('fraction of cue cells');
xlabel('correlation with value model');
yticks(0:0.1:0.3);

%trial type
tt=valueShuffles(includedModels,1)==valueShuffles(includedModels,2) &...
    valueShuffles(includedModels,3)==valueShuffles(includedModels,4) &...
    valueShuffles(includedModels,5)==valueShuffles(includedModels,6);
subplot(4,5,7);
hold on
scatter(corrWVal(tt==0,1),frequency(tt==0),24,[0 0 0],'filled');
scatter(corrWVal(tt,1),frequency(tt),24,[1 0 0.3],'filled');
plot([-1 1],[1/length(includedModels) 1/length(includedModels)],':','color','k','linewidth',0.5);
ylim([0 0.2]);
ylabel('fraction of cue cells');
xlabel('correlation with value model');
yticks(0:0.1:0.3);


category=9*ones(length(improvement),1);
category(cueK&~behavior)=4; %cue neurons, other odor coding schemes
category(cueK&~behavior&bestModel==1)=1; %value
category(cueK&~behavior&ismember(bestModel,valRank(2:17)))=2; %value like
category(cueK&~behavior&bestModel==max(bestModel))=3; %untuned




subplot(4,5,8);
hold on
valcand=ismember(bestModel,[1 2]);
scatter(varExpValueNew(valcand,1),varExpValueNew(valcand,91));
plot([0 1],[0 1],'color','k');
xlabel('variance explained with 1,0.7,0.1');
ylabel('variance explained with 1,0.5,0');
title('value cells (selected with either)');

% visual depiction of shuffles
colormap('summer');
subplot(5,1,4);
imagesc(1:length(includedModels),1:6,valueShuffles(91:end,:)',[0 1]);
xticks([]);
yticks(1:6);
yticklabels({'CS+','CS+','CS50','CS50','CS-','CS-'});
title('Value assigned to each odor in each shuffle');
colorbar;

% visual depiction of shuffles, reordered by correlation
subplot(5,1,5);
includedShuffles = valueShuffles(91:end,:)';

imagesc(1:length(includedModels),1:6,includedShuffles(:,valRank),[0 1]);
xticks([]);
yticks(1:6);
yticklabels({'CS+','CS+','CS50','CS50','CS-','CS-'});
title('Value assigned to each odor in each shuffle, ordered by similarity to value');
colorbar;

%% get proportions of value and meaning and compare between regions
regions={'ALM','ACA','FRP','PL','ILA','ORB','DP','TTd','AON'};
reggrps=[2:4];
neuronGroupAdj = neuronGroup + 1;
neuronGroupAdj(neuronGroupAdj==5)=1;
pcutoff=0.05/36;

% regions={'ACA','FRP','PL','ILA','ORB','CP','ACB'};
% reggrps=[2 4];
neuronRegionOlf=neuronRegionAdj;
neuronRegionOlf(ismember(neuronRegionAdj,{'EPd','PIR'}))={'OLF'};
allReg=unique(neuronRegionOlf);
order=[9 10 5 2 3 1 7 12 8 11 6 13 4];
allReg=allReg(order);

regInd=NaN(length(regions),1);
for reg=1:length(regions)
    regInd(reg,1)=find(ismember(allReg,regions(reg)));
end


%fit all together, allowing random effect of session to have more power
catcolors={[0 0.7 1],[0.7 0 1],[0.6 0.6 0.6],[0.05 0.25 0.45]};
titles={'value','value-like','untuned','other'};
cats={1,2,3,4};
for ct=1:length(cats)
    %region
    barData=zeros(2,5);
    subplot(4,4,ct);
    hold on;
    sel=ismember(category,cats{ct});
    othersel=ismember(category,1:8)&~sel;
    regionnumber=zeros(totalNeurons,1);
    for reg=1:length(allReg)
        regsel=ismember(neuronRegionOlf,allReg(reg));
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel);
        barData(2,reg)=sum(othersel&regsel)/sum(regsel);   
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial','dummyvarcoding','reference');
    
    [ypred,ypredCI,df]=predict(lm,'conditional',false);
        
    Hall=[];
    Htemp = [0 0 0 0 0 0 0 0 0 0 0 0 0];
    for reg1 = 1:length(regions)-1
        for reg2 = reg1+1:length(regions)
            H=Htemp;
            H(reg1+4)=1;
            H(reg2+4)=-1;
            Hall=cat(1,Hall,H);
        end
    end
    greater=zeros(length(regions));
    for reg1 = 1:length(regions)
        for reg2 = 1:length(regions)
            H=Htemp;
            H(reg1+4)=1;
            H(reg2+4)=-1;
            greater(reg1,reg2)=coefTest(lm,H);
        end
    end    
    
    %Hall=unique(Hall,'rows');
    pVals=NaN(length(Hall),1);
    for p=1:length(pVals)
        pVals(p,1) = coefTest(lm,Hall(p,:));
    end
    barData=barData(:,regInd);
    b=bar(barData','stacked');
    %b(3).FaceColor=[0.4 0 0.5]; %odor
    b(2).FaceColor=[0.85 0.85 0.85]; %meaning
    b(1).FaceColor=catcolors{ct}; %value
    %b(4).FaceColor=[0 0 0.2]; %odor
    % b(5).FaceColor=[0.9 0.2 0];
    %b(3).FaceColor=[1 1 1];
    upperCI=[];
    lowerCI=[];
    estmean=[];
    for reg=1:length(regions)
        regsel=find(ismember(neuronRegionOlf,regions(reg)),1);
        estmean(1,reg)=ypred(regsel);
        upperCI(1,reg)=ypredCI(regsel,2);
        lowerCI(1,reg)=ypredCI(regsel,1);
    end
    
    errorbar(1:length(regions),estmean(1,:),estmean(1,:)-lowerCI(1,:),upperCI(1,:)-estmean(1,:),'o','linewidth',0.75,'color','k');
    
    xlim([0.5 length(regions)+0.5]);
    %legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
    xticks(1:length(regions));
    xticklabels(regions);
    xtickangle(45);
    ylim([0 0.55]);
    
    
    %make heatmap
    target=catcolors{ct};
    customMap=NaN(100,3);
    for i=1:length(customMap)
        customMap(i,:)=target*i/(length(customMap));
    end
    customMap=cat(1,[1 1 1],customMap);
    
    %by region
    %greater=zeros(length(regions));
    frac=zeros(length(regions));
    for r1=1:length(regions)
        for r2=1:length(regions)
            frac(r1,r2)=barData(1,r1)-barData(1,r2);
            %greater(r1,r2)=lowerCI(1,r1)>upperCI(1,r2);
        end        
    end
    s=subplot(4,4,4+ct);
    imagesc(frac .* (greater<pcutoff & frac>0), [0 ceil(max(frac,[],'all')*20)/20]);
    
    colormap(s,customMap);
    cb=colorbar;
    set(cb,'ticks',[0:0.1:0.3]);
    xticks(1:length(regions));
    xticklabels(regions);
    xtickangle(45);
    yticks(1:length(regions));
    yticklabels(regions);
    title(titles{ct});
    if ct==1
        ylabel('increase in this region...');
        xlabel('over this region');
    end    
    
    %region group
    barData=zeros(2,3);
    subplot(4,4,8+ct);
    hold on;
    sel=ismember(category,cats{ct});
    othersel=ismember(category,1:8)&~sel;    
    regionnumber=zeros(totalNeurons,1);
    for reg=1:4
        regsel=ismember(neuronGroupAdj,reg);
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel);   
        barData(2,reg)=sum(othersel&regsel)/sum(regsel);   
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial','dummyvarcoding','reference');
    [ypred,ypredCI]=predict(lm,'conditional',false);
    
    Htemp = [0 0 0 0];
    greater=zeros(length(reggrps));
    for reg1 = 1:length(reggrps)
        for reg2 = 1:length(reggrps)
            H=Htemp;
            H(reg1+1)=1;
            H(reg2+1)=-1;
            greater(reg1,reg2)=coefTest(lm,H);
        end
    end    
    
    barData=barData(:,reggrps);    
    b=bar(barData','stacked');
    %b(3).FaceColor=[0.4 0 0.5]; %odor
    b(2).FaceColor=[0.85 0.85 0.85];
    b(1).FaceColor=catcolors{ct}; %value
    %b(4).FaceColor=[0 0 0.2]; %odor
    % b(5).FaceColor=[0.9 0.2 0];
    %b(3).FaceColor=[1 1 1];
    upperCI=[];
    lowerCI=[];
    estmean=[];
    for reg=1:length(reggrps)
        regsel=find(ismember(neuronGroupAdj,reggrps(reg)),1);
        estmean(1,reg)=ypred(regsel);
        upperCI(1,reg)=ypredCI(regsel,2);
        lowerCI(1,reg)=ypredCI(regsel,1);
    end
    
    errorbar(1:length(reggrps),estmean(1,:),estmean(1,:)-lowerCI(1,:),upperCI(1,:)-estmean(1,:),'o','linewidth',0.75,'color','k');
    
    
    Htemp = [0 0 0 0];
    greater=zeros(length(reggrps));
    for reg1 = 1:length(reggrps)
        for reg2 = 1:length(reggrps)
            H=Htemp;
            H(reg1+1)=1;
            H(reg2+1)=-1;
            greater(reg1,reg2)=coefTest(lm,H);
        end
    end      
    
    
    xlim([0.5 7+0.5]);
    %legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
    xticks(1:length(reggrps));
    xticklabels(regionGroupNames(reggrps));
    xtickangle(45);
    ylim([0 0.55]);
    
    
    %make heatmap
    target=catcolors{ct};
    customMap=NaN(100,3);
    for i=1:length(customMap)
        customMap(i,:)=target*i/(length(customMap));
    end
    customMap=cat(1,[1 1 1],customMap);
    
    %by region
    %greater=zeros(length(reggrps));
    frac=zeros(length(reggrps));
    for r1=1:length(reggrps)
        for r2=1:length(reggrps)
            frac(r1,r2)=barData(1,r1)-barData(1,r2);
            %greater(r1,r2)=lowerCI(1,r1)>upperCI(1,r2);
        end        
    end
    s=subplot(4,4,12+ct);
    imagesc(frac .* (greater<0.05/3 & frac>0), [0 ceil(max(frac,[],'all')*20)/20]);
    
    colormap(s,customMap);
    cb=colorbar;
    set(cb,'ticks',[0:0.1:0.3]);
    xticks(1:length(reggrps));
    xticklabels(regionGroupNames(reggrps-1));
    xtickangle(45);
    yticks(1:length(reggrps));
    yticklabels(regionGroupNames(reggrps-1));
    title(titles{ct});
    if ct==1
        ylabel('increase in this region...');
        xlabel('over this region');
    end        
    
end


%% get proportions of value and meaning out of cue cells
regions={'ALM','ACA','FRP','PL','ILA','ORB','DP','TTd','AON'};
reggrps=[2:4];
neuronGroupAdj = neuronGroup + 1;
neuronGroupAdj(neuronGroupAdj==5)=1;
pcutoff=0.05/36;

%regions={'ACA','FRP','PL','ILA','ORB','CP','ACB'};
%reggrps=[2 4];

allReg=unique(neuronRegionOlf);
order=[9 10 5 2 3 1 7 12 8 11 6 13 4];
allReg=allReg(order);
regInd=NaN(length(regions),1);
for reg=1:length(regions)
    regInd(reg,1)=find(ismember(allReg,regions(reg)));
end

%fit all together, allowing random effect of session to have more power
catcolors={[0 0.7 1],[0.7 0 1],[0.6 0.6 0.6],[0.05 0.25 0.45]};
titles={'value','value-like','untuned','other'};
cats={1,2,3,4};
for ct=1:length(cats)
    %region
    barData=zeros(2,5);
    subplot(4,4,ct);
    hold on;
    sel=ismember(category,cats{ct});
    othersel=ismember(category,1:8)&~sel;
    regionnumber=zeros(totalNeurons,1);
    for reg=1:length(allReg)
        regsel=ismember(neuronRegionOlf,allReg(reg))&category<9;
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel);
        barData(2,reg)=sum(othersel&regsel)/sum(regsel);   
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial');
    [ypred,ypredCI]=predict(lm,'conditional',false);
    
    
    Hall=[];
    Htemp = [0 0 0 0 0 0 0 0 0 0 0 0 0];
    for reg1 = 1:length(regions)-1
        for reg2 = reg1+1:length(regions)
            H=Htemp;
            H(reg1+4)=1;
            H(reg2+4)=-1;
            Hall=cat(1,Hall,H);
        end
    end
    greater=zeros(length(regions));
    for reg1 = 1:length(regions)
        for reg2 = 1:length(regions)
            H=Htemp;
            H(reg1+4)=1;
            H(reg2+4)=-1;
            greater(reg1,reg2)=coefTest(lm,H);
        end
    end    
    
    %Hall=unique(Hall,'rows');
    pVals=NaN(length(Hall),1);
    for p=1:length(pVals)
        pVals(p,1) = coefTest(lm,Hall(p,:));
    end    
    
    
    barData=barData(:,regInd);
    b=bar(barData','stacked');
    %b(3).FaceColor=[0.4 0 0.5]; %odor
    b(2).FaceColor=[0.85 0.85 0.85]; %meaning
    b(1).FaceColor=catcolors{ct}; %value
    %b(4).FaceColor=[0 0 0.2]; %odor
    % b(5).FaceColor=[0.9 0.2 0];
    %b(3).FaceColor=[1 1 1];
    upperCI=[];
    lowerCI=[];
    estmean=[];
    for reg=1:length(regions)
        regsel=find(ismember(neuronRegionOlf(category<9),regions(reg)),1);
        estmean(1,reg)=ypred(regsel);
        upperCI(1,reg)=ypredCI(regsel,2);
        lowerCI(1,reg)=ypredCI(regsel,1);
    end
     
    errorbar(1:length(regions),estmean(1,:),estmean(1,:)-lowerCI(1,:),upperCI(1,:)-estmean(1,:),'o','linewidth',0.75,'color','k');
    
    xlim([0.5 length(regions)+0.5]);
    %legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
    xticks(1:length(regions));
    xticklabels(regions);
    xtickangle(45);
    ylim([0 1]);
    
    
    %make heatmap
    target=catcolors{ct};
    customMap=NaN(100,3);
    for i=1:length(customMap)
        customMap(i,:)=target*i/(length(customMap));
    end
    customMap=cat(1,[1 1 1],customMap);
    
    %by region
    %greater=zeros(length(regions));
    frac=zeros(length(regions));
    for r1=1:length(regions)
        for r2=1:length(regions)
            frac(r1,r2)=barData(1,r1)-barData(1,r2);
            %greater(r1,r2)=lowerCI(1,r1)>upperCI(1,r2);
        end        
    end
    s=subplot(4,4,4+ct);
    imagesc(frac .* (greater<pcutoff & frac>0), [0 ceil(max(frac,[],'all')*20)/20]);
    
    colormap(s,customMap);
    cb=colorbar;
    set(cb,'ticks',[0:0.1:0.3]);
    xticks(1:length(regions));
    xticklabels(regions);
    xtickangle(45);
    yticks(1:length(regions));
    yticklabels(regions);
    title(titles{ct});
    if ct==1
        ylabel('increase in this region...');
        xlabel('over this region');
    end    
 
    
    %region group
    barData=zeros(2,3);
    subplot(4,4,8+ct);
    hold on;
    sel=ismember(category,cats{ct});
    othersel=ismember(category,1:8)&~sel;    
    regionnumber=zeros(totalNeurons,1);
    for reg=1:4
        regsel=ismember(neuronGroupAdj,reg)&category<9;
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel); 
        barData(2,reg)=sum(othersel&regsel)/sum(regsel);   
        
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial');
    [ypred,ypredCI]=predict(lm,'conditional',false);

   
    Htemp = [0 0 0 0];
    greater=zeros(length(reggrps));
    for reg1 = 1:length(reggrps)
        for reg2 = 1:length(reggrps)
            H=Htemp;
            H(reg1+1)=1;
            H(reg2+1)=-1;
            greater(reg1,reg2)=coefTest(lm,H);
        end
    end     
     
    
    barData=barData(:,reggrps);    
    b=bar(barData','stacked');
    %b(3).FaceColor=[0.4 0 0.5]; %odor
    b(2).FaceColor=[0.85 0.85 0.85];
    b(1).FaceColor=catcolors{ct}; %value
    %b(4).FaceColor=[0 0 0.2]; %odor
    % b(5).FaceColor=[0.9 0.2 0];
    %b(3).FaceColor=[1 1 1];
    upperCI=[];
    lowerCI=[];
    estmean=[];
    for reg=1:length(reggrps)
        regsel=find(ismember(neuronGroupAdj(category<9),reggrps(reg)),1);
        estmean(1,reg)=ypred(regsel);
        upperCI(1,reg)=ypredCI(regsel,2);
        lowerCI(1,reg)=ypredCI(regsel,1);
    end
    
    errorbar(1:length(reggrps),estmean(1,:),estmean(1,:)-lowerCI(1,:),upperCI(1,:)-estmean(1,:),'o','linewidth',0.75,'color','k');
    
    xlim([0.5 7+0.5]);
    %legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
    xticks(1:length(reggrps));
    xticklabels(regionGroupNames(reggrps-1));
    xtickangle(45);
    ylim([0 1]);
    
    
    %make heatmap
    target=catcolors{ct};
    customMap=NaN(100,3);
    for i=1:length(customMap)
        customMap(i,:)=target*i/(length(customMap));
    end
    customMap=cat(1,[1 1 1],customMap);
    
    %by region
    %greater=zeros(length(reggrps));
    frac=zeros(length(reggrps));
    for r1=1:length(reggrps)
        for r2=1:length(reggrps)
            frac(r1,r2)=barData(1,r1)-barData(1,r2);
            %greater(r1,r2)=lowerCI(1,r1)>upperCI(1,r2);
        end        
    end
    s=subplot(4,4,12+ct);
    imagesc(frac .* (greater<0.05/3 & frac>0), [0 ceil(max(frac,[],'all')*20)/20]);
    
    colormap(s,customMap);
    cb=colorbar;
    set(cb,'ticks',[0:0.1:0.3]);
    xticks(1:length(reggrps));
    xticklabels(regionGroupNames(reggrps-1));
    xtickangle(45);
    yticks(1:length(reggrps));
    yticklabels(regionGroupNames(reggrps-1));
    title(titles{ct});
    if ct==1
        ylabel('increase in this region...');
        xlabel('over this region');
    end        
    
end


%% value models schematic
examplePerms=[91 123];
exampleNeurons=[3 10 20];

plotWindow=[0 2.5];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);

binTimesGLMTrl = analysisWindow(1)+binSize/2:binSize:analysisWindow(2)-binSize/2;
plotBinsGLM=binTimesGLMTrl>=plotWindow(1) & binTimesGLMTrl<=plotWindow(2);


colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];
specs={'-','--'};

folds=4;

neuronNumber=[1:length(neuronSession)]';
%get neuron number from order in value category
entries=find(category==1);
for n=1:length(exampleNeurons)
    nn=entries(exampleNeurons(n));
    ta=[];
    subplot(length(exampleNeurons),5,(n-1)*5+1);
    %subplot(6,6,n);
    
    for os=1:2
        
        
        hold on;
        for cue=1:3
            activity=PSTH3{cue,os}(nn,plotBins);
            ta=[ta activity];
            plot(binTimes(plotBins),activity,specs{os},'linewidth',1,'color',colors{cue});
            
        end
        
    end
    plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
    axis([plotWindow floor(min(ta)) ceil(max(ta))])
    if n>1
        xticks([]);
    end
    if n==1
        xticks([0 2.5]);
        xlabel('seconds from odor onset');
        ylabel('z-score');
    end
    title(num2str(exampleNeurons(n)));
    
    %get model fits
    session=neuronSession(nn);
    sessionNeuron=find(neuronNumber(neuronSession==session)==nn);
    sessionFolder=fullfile(direc,sessionSubject{session},sessionDate{session});
    cueTimes = readNPY(fullfile(sessionFolder,'cues.times.npy'));
    rewprob = readNPY(fullfile(sessionFolder,'cues.rewardProbability.npy'));
    odorset = readNPY(fullfile(sessionFolder,'cues.odorSet.npy'));
    cueOrder=zeros(length(cueTimes),1);
    cueOrder(rewprob==1&odorset==1)=1;
    cueOrder(rewprob==1&odorset==2)=2;
    cueOrder(rewprob==0.5&odorset==1)=3;
    cueOrder(rewprob==0.5&odorset==2)=4;
    cueOrder(rewprob==0&odorset==1)=5;
    cueOrder(rewprob==0&odorset==2)=6;
    
    A = Xall{session}(:,[1:binsPerCue size(Xall{session},2)-5:size(Xall{session},2)]);
    
    for s=1:length(examplePerms)
        vs=examplePerms(s);
        for c=1:length(cueTimes)
            cueValue=valueShuffles(vs,cueOrder(c));
            dm=cueValue*eye(binsPerCue);
            A(c*binsPerTrial-binsPerCue+1:c*binsPerTrial,1:binsPerCue)=dm;
        end
        y=Yall{session}(:,sessionNeuron);
        fit = glmnet(A,y,[],glmopts);
        kernel = fit.beta(1:sum(plotBinsGLM),end);
        
        
        
        subplot(length(exampleNeurons),5,(n-1)*5+1+s);
        hold on;
        o=0;
        for cue=1:3
            for os=1:2
                o=o+1;
                plot(binTimesGLMTrl(plotBinsGLM),kernel(1:sum(plotBinsGLM))*valueShuffles(vs,o),specs{os},'linewidth',1,'color',colors{cue});
            end           
        end
        plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
        axis([plotWindow floor(min(ta)) ceil(max(ta))])
        if n>1
            xticks([]);
        end
        if n==1
            title(sprintf('shuffle %d',examplePerms(s)));
        end
        
    end
    
end


%% activity of cue cell types
figure;
map=redblue(256);

cueWindow=[0 2.5];
cueBins=binTimes>=cueWindow(1) & binTimes<=cueWindow(2);

%heatmaps
plotWindow=[-0.5 2.5];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);

sortWindow=[0 2.5];
sortBins=binTimes>=sortWindow(1) & binTimes<=sortWindow(2);

colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];

pn=0;
activity=[];
thisActivity=PSTH3{1,1};
cueResp=mean(thisActivity(:,sortBins),2);
cats={1,2,3};

% regions={'ALM','MOs','ACA','FRP','PL','ILA','ORB','CP','ACB','DP','TTd','AON','OLF'};
% 
% neuronRegionOlf=neuronRegionAdj;
% neuronRegionOlf(ismember(neuronRegionAdj,{'EPd','PIR'}))={'OLF'};
% regionCode=NaN(length(neuronRegionOlf),1);
% for NN=1:length(neuronRegionOlf)
%     regionCode(NN,1)=find(ismember(regions,neuronRegionOlf(NN)));
% end
% regionCodeOpp=12-regionCode;
regionCodeOpp=ones(totalNeurons,1);

for ct=1:length(cats)
sel=ismember(category,cats{ct});
sortcrit=[regionCodeOpp(sel) cueResp(sel)];
[~,sortOrder]=sortrows(sortcrit,[1 2]);
pn=0;
for cue=1:3
    for os=1:2
        
        pn=pn+1;
        
        ax(1)=subplot(3,15,(ct-1)*15+pn);
        colormap(ax(1),map);
        hold on;
        
        
        activity = PSTH3{cue,os}(sel,:);
        activity = activity(sortOrder,:);
        

        
        imagesc(binTimes(plotBins),[1 length(activity)],activity(:,plotBins),[-3 3]);%[min(min(cspActivity)) max(max(cspActivity))]);
        ylim([0.5 size(activity,1)+0.5])
        plot([0 0],[0.5 sum(sel)+0.5],'color',colors{cue},'linewidth',0.75);
        set(gca,'ytick',[]);
        if pn==1 ylabel(sprintf('n=%d',sum(sel))); end
        if pn==1 xlabel('seconds from cue'); end
        if pn>1 xticks([]); end
        

        
    end
    
    
    
end

colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];
directions=[1 -1];
specs={'-','--'};
diractivity = mean(cat(3,PSTH6{(1-1)*2+2,1},PSTH6{(1-1)*2+1,2}),3);
direction=sign(max(diractivity,[],2)-abs(min(diractivity,[],2)));

    for d=1:2
        subplot(6,3,3+(ct-1)*6+(d-1)*3);
        
        hold on;
        for cue=1:3
            for os=1:2
                
                activity = PSTH6{(cue-1)*2+os,os};
                
                %get values
                psth=nanmean(activity(sel&direction==directions(d),plotBins));
                sem=nanste(activity(sel&direction==directions(d),plotBins),1); %calculate standard error of the mean
                up=psth+sem;
                down=psth-sem;
                
                %plotting
                plot(binTimes(plotBins),psth,specs{os},'Color',colors{cue,1},'linewidth',0.75);
                %patch([xvals,xvals(end:-1:1)],[up,down(end:-1:1)],regcolors{reg,1},'EdgeColor','none');alpha(0.2);
                
                
            end
            
        end
        
        if d==1 text(0.1,1.5,num2str(sum(sel&direction==directions(d)))); end
        if d==2 text(0.1,-0.45,num2str(sum(sel&direction==directions(d)))); end
        if reg==1 & d==2
            xlabel('seconds from odor onset');
            ylabel('z-score');
        else
            xticks([]);
        end
        
        if reg>1 yticks([]); end
        plot([0 0],[-1 2.5],'color','k','linewidth',0.75);
        plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
        if d==1 axis([plotWindow -0.25 2.5]); end
        if d==2 axis([plotWindow -1 0.2]); end
    end




end

%%  plot PC1 value from each region
figure;

cueWindow=[0 2.5];
cueBins=binTimes>=cueWindow(1) & binTimes<=cueWindow(2);

neuronRegionOlf=neuronRegionAdj;
neuronRegionOlf(ismember(neuronRegionAdj,{'EPd','PIR'}))={'OLF'};
division=neuronRegionOlf;

regions={'ALM','ACA','FRP','PL','ILA','ORB','DP','TTd','AON','CP','ACB'};

plotWindow=[-0.5 2.5];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);
for reg=1:length(regions)
    sel=ismember(division,regions(reg)) & ismember(category,[2]);
activity=[];
for cue=1:3
    for os=1:2   
        activity = [activity PSTH3{cue,os}(sel,plotBins)];
    end
end
activity=activity./max(abs(activity),[],2);
[coeff,score,~,~,explained]=pca(activity');
specs={'-','--'};
comp=1;
if strcmp(regions(reg),'ACA') comp=1; end
cond=0;
if explained(comp)>5
    subplot(4,9,reg)
    
    hold on
    for cue=1:3
        for os=1:2
            cond=cond+1;
            plot(binTimes(plotBins),-score((cond-1)*sum(plotBins)+1:cond*sum(plotBins),comp),specs{os},'color',colors{cue},'linewidth',1);
        end
    end
    ylabel(sprintf('#%g (%g%%)',comp,round(explained(comp),1)));
    yticks([]);
    xticks([0 2.5]);
    title(regions{reg});
    
    %     plot([0 0],[min(score(:,comp)) max(score(:,comp))],':','color','k','linewidth',0.5);
    %     plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
    
end
xlim(cueWindow);
xlabel(sprintf('n = %d',sum(sel)));



end

division=neuronGroup;
for reg=1:4
    sel=ismember(division,reg) & category==2;
activity=[];
for cue=1:3
    for os=1:2   
        activity = [activity PSTH3{cue,os}(sel,plotBins)];
    end
end
activity=activity./max(abs(activity),[],2);
[coeff,score,~,~,explained]=pca(activity');
specs={'-','--'};
for comp=1
    cond=0;
    if explained(comp)>5
    subplot(4,9,18+reg)
    
    hold on
    for cue=1:3
        for os=1:2
            cond=cond+1;
            plot(binTimes(plotBins),-score((cond-1)*sum(plotBins)+1:cond*sum(plotBins),comp),specs{os},'color',colors{cue},'linewidth',1);
        end
    end
    ylabel(sprintf('#%g (%g%%)',comp,round(explained(comp),1)));
    yticks([]);
    xticks([0 2.5]);
    xlabel(sprintf('n = %d',sum(sel)));
    if comp==1 title(regionGroupNames{reg}); end
    
%     plot([0 0],[min(score(:,comp)) max(score(:,comp))],':','color','k','linewidth',0.5);
%     plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
    
    end
    xlim(cueWindow);
end

end



%% GLM spatial plotting, 3D
figure;
catcolors={[0 0.7 1],[0.7 0 1],[0.6 0.6 0.6]};

%all
neuronXYZmm=neuronXYZ/1000;
cats={1,2,3};
for ct=1:3
    subplot(1,3,ct);
    hold on;
    sel=ismember(category,cats{ct});
    offset=(rand(sum(sel),1)*10-5)/1000;
    p1=scatter3(abs(neuronXYZmm(sel,1))+offset,neuronXYZmm(sel,2)+offset,neuronXYZmm(sel,3)+offset,24,catcolors{ct},'filled');
    p1.MarkerFaceAlpha=0.2;
view(20,20);
axis([0 3 1 3 -7 -1]);
if ct==1
xlabel('ML (mm)');
ylabel('AP (mm)');
zlabel('DV (mm)');
else
    xticks([]);
    yticks([]);
    zticks([]);
end

end


%% coding dimension for CS+/CS50/CS-
figure;
%categories for subselection
sels={};
sels{1}=ismember(category,3);
sels{2}=ismember(category,2);
sels{3}=category==1;
colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];
vcolors={[0.6 0.6 0.6],[0.7 0 1],[0 0.7 1]};

numstraps=5000;
amin=-0.1;
amax=1.5;
cueWindow=[0 2.5];
cueBins=binTimes>=cueWindow(1) & binTimes<=cueWindow(2);

plotWindow=[-1 10];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);
plotBinTimes=binTimes(plotBins);
baseBins=plotBinTimes<0;

%normalize activity from 0-1
dffPSTH6n={};
dffPSTH3n={};
for session=1
    dffPSTH6all=[];
    for condition=1:6
        dffPSTH6all=cat(2,dffPSTH6all,PSTH6{condition,session}(:,cueBins));
    end
end
for session=1:2
    for condition=1:6
       dffPSTH6n{condition,session}=PSTH6{condition,session}./max(abs(dffPSTH6all),[],2); 
    end
end

%get activity for dimensions
csmActivity=dffPSTH6n{5,1}(:,cueBins);
cspActivity=dffPSTH6n{1,1}(:,cueBins);
csfActivity=dffPSTH6n{3,1}(:,cueBins);

bigBinSize=5;
bigBins=floor(sum(cueBins)/bigBinSize);
bsf=0;
for bb=1:bigBins
    csmActivityB(:,bb)=mean(csmActivity(:,bsf+1:bsf+bigBinSize),2);
    cspActivityB(:,bb)=mean(cspActivity(:,bsf+1:bsf+bigBinSize),2);
    csfActivityB(:,bb)=mean(csfActivity(:,bsf+1:bsf+bigBinSize),2);
    
    bsf=bsf+bigBinSize;
end

[~,PMind]=max(abs(cspActivityB-csmActivityB),[],2,'linear');
[~,PFind]=max(abs(cspActivityB-csfActivityB),[],2,'linear');
[~,FMind]=max(abs(csfActivityB-csmActivityB),[],2,'linear');

differenceVectorPM=cspActivityB(PMind)-csmActivityB(PMind);
differenceVectorPF=cspActivityB(PFind)-csfActivityB(PFind);
differenceVectorFM=csfActivityB(FMind)-csmActivityB(FMind);


sproj={};
oproj={};
disttime={};

%for getting 0 and 1, use z-score
csmActivity=PSTH6{5,1}(:,cueBins);
cspActivity=PSTH6{1,1}(:,cueBins);
csfActivity=PSTH6{3,1}(:,cueBins);

bsf=0;
for bb=1:bigBins
    csmActivityB(:,bb)=mean(csmActivity(:,bsf+1:bsf+bigBinSize),2);
    cspActivityB(:,bb)=mean(cspActivity(:,bsf+1:bsf+bigBinSize),2);
    csfActivityB(:,bb)=mean(csfActivity(:,bsf+1:bsf+bigBinSize),2);
    
    bsf=bsf+bigBinSize;
end

projCorrStrp=NaN(numstraps,length(sels));
baseCSMdist=NaN(numstraps,length(sels));
cspfangle=NaN(numstraps,length(sels));
for cl=1:length(sels)
    sel=sels{cl};
    
    PF1=differenceVectorPF(sel)' * cspActivityB(PFind(sel));
    PF0=differenceVectorPF(sel)' * csfActivityB(PFind(sel));
    PM1=differenceVectorPM(sel)' * cspActivityB(PMind(sel));
    PM0=differenceVectorPM(sel)' * csmActivityB(PMind(sel));
    FM1=differenceVectorFM(sel)' * csfActivityB(FMind(sel));
    FM0=differenceVectorFM(sel)' * csmActivityB(FMind(sel));

    for cue=1:3
        sproj{cue,1}=(differenceVectorPM(sel)' * PSTH6{(cue-1)*2+2,1}(sel,:) - PM0) ./ (PM1 - PM0);
        sproj{cue,2}=(differenceVectorPF(sel)' * PSTH6{(cue-1)*2+2,1}(sel,:) - PF0) ./ (PF1 - PF0);
        sproj{cue,3}=(differenceVectorFM(sel)' * PSTH6{(cue-1)*2+2,1}(sel,:) - FM0) ./ (FM1 - FM0);
        
        oproj{cue,1}=(differenceVectorPM(sel)' * PSTH3{cue,2}(sel,:) - PM0) ./ (PM1 - PM0);
        oproj{cue,2}=(differenceVectorPF(sel)' * PSTH3{cue,2}(sel,:) - PF0) ./ (PF1 - PF0);
        oproj{cue,3}=(differenceVectorFM(sel)' * PSTH3{cue,2}(sel,:) - FM0) ./ (FM1 - FM0);
    end
    
    subplot(3,4,1+(cl-1)*4);
    hold on;
    for cue=1:3
        plot(sproj{cue,1}(plotBins),sproj{cue,3}(plotBins),':','linewidth',1,'color',colors{cue});
        plot(sproj{cue,1}(cueBins),sproj{cue,3}(cueBins),'linewidth',1,'color',colors{cue});
    end
    axis([amin amax amin amax amin amax]);
    if cl==1
        title('same odor set');
            xlabel('CS- ---> CS+');
        ylabel('CS- ---> CS50');
    end
    xticks([0 1]);
    yticks([0 1]);
    
    subplot(3,4,2+(cl-1)*4);
    hold on;
    for cue=1:3
%        plot(oproj{cue,1}(plotBins),oproj{cue,3}(plotBins),':','linewidth',1,'color',colors{cue});
        plot(oproj{cue,1}(cueBins),oproj{cue,3}(cueBins),'linewidth',1,'color',colors{cue});
    
    end
    axis([amin amax amin amax amin amax]);
    xticks([]);
    yticks([]);
    if cl==1 title('other odor set'); end     
     
     
     %bootstrap a statistic
     incNeur=find(sel);

     for strap=1:numstraps
         ns=incNeur(randsample(length(incNeur),length(incNeur),'true'));
         PF1=differenceVectorPF(ns)' * cspActivityB(PFind(ns));
         PF0=differenceVectorPF(ns)' * csfActivityB(PFind(ns));
         PM1=differenceVectorPM(ns)' * cspActivityB(PMind(ns));
         PM0=differenceVectorPM(ns)' * csmActivityB(PMind(ns));
         FM1=differenceVectorFM(ns)' * csfActivityB(FMind(ns));
         FM0=differenceVectorFM(ns)' * csmActivityB(FMind(ns));
         
         for cue=1:3
             sproj{cue,1}=(differenceVectorPM(ns)' * PSTH6{(cue-1)*2+2,1}(ns,cueBins) - PM0) ./ (PM1 - PM0);
             sproj{cue,2}=(differenceVectorPF(ns)' * PSTH6{(cue-1)*2+2,1}(ns,cueBins) - PF0) ./ (PF1 - PF0);
             sproj{cue,3}=(differenceVectorFM(ns)' * PSTH6{(cue-1)*2+2,1}(ns,cueBins) - FM0) ./ (FM1 - FM0);
             
             oproj{cue,1}=(differenceVectorPM(ns)' * PSTH3{cue,2}(ns,cueBins) - PM0) ./ (PM1 - PM0);
             oproj{cue,2}=(differenceVectorPF(ns)' * PSTH3{cue,2}(ns,cueBins) - PF0) ./ (PF1 - PF0);
             oproj{cue,3}=(differenceVectorFM(ns)' * PSTH3{cue,2}(ns,cueBins) - FM0) ./ (FM1 - FM0);
         end
     
         %only use CS- subtracted axes for correlation
         sprojall=cat(1,sproj{:,[1 3]});
         oprojall=cat(1,oproj{:,[1 3]});
         projCorrStrp(strap,cl)=corr(sprojall(:),oprojall(:));
         
         %baseline distance from CS-
         basex=mean([sproj{1,1}(baseBins) sproj{2,1}(baseBins) sproj{3,1}(baseBins)]);
         basey=mean([sproj{1,3}(baseBins) sproj{2,3}(baseBins) sproj{3,3}(baseBins)]);
         baseCSMdist(strap,cl)=sqrt(basex^2+basey^2);
         
         %angle between CS+ at peak along (CS-/CS+) and CS50 at peak along (CS-/CS50)
         [cspx,bin]=max(sproj{1,1}); %CSP
         cspy=sproj{1,3}(bin);        
         [csfy,bin]=max(sproj{2,3}); %CSF
         csfx=sproj{2,1}(bin);        
         v_1 = [cspx,cspy,0] - [basex,basey,0];
         v_2 = [csfx,csfy,0] - [basex,basey,0];
         cspfangle(strap,cl) = atan2(norm(cross(v_1, v_2)), dot(v_1, v_2)) * 180/pi; %angle
         
     end
     

     

end

subplot(2,4,7);
hold on;
for c=1:3
histogram(baseCSMdist(:,c),0.1:0.01:1,'facecolor',vcolors{c},'edgecolor',vcolors{c},'orientation','horizontal','normalization','probability');
end

ylabel('dist from cs-');
xlabel('frequency');
legend('odor','meaning','value');
bootp=(sum(sum(baseCSMdist(:,1)<=baseCSMdist(:,3)'))+1)/(numstraps^2+1);
text(0.1,0.9,sprintf('p = %g',round(bootp,2,'significant')));
bootp=(sum(sum(baseCSMdist(:,2)<=baseCSMdist(:,3)'))+1)/(numstraps^2+1);
text(0.1,0.8,sprintf('p = %g',round(bootp,2,'significant')));
ylim([0 0.8]);

subplot(2,4,8);
hold on;
for c=1:3
histogram(cspfangle(:,c),0:1:110,'facecolor',vcolors{c},'edgecolor',vcolors{c},'orientation','horizontal','normalization','probability');
end

ylabel('CS+, CS50 angle');
xlabel('frequency');
legend('odor','meaning','value');
bootp=(sum(sum(cspfangle(:,1)<=cspfangle(:,3)'))+1)/(numstraps^2+1);
text(0.05,25,sprintf('p = %g',round(bootp,2,'significant')));
bootp=(sum(sum(cspfangle(:,2)<=cspfangle(:,3)'))+1)/(numstraps^2+1);
text(0.15,25,sprintf('p = %g',round(bootp,2,'significant')));
bootp=(sum(sum(cspfangle(:,1)>=cspfangle(:,2)'))+1)/(numstraps^2+1);
text(0.05,95,sprintf('p = %g',round(bootp,2,'significant')));
ylim([0 45]);
yticks([0:15:45]);


subplot(2,4,4);
hold on;
for c=1:3
histogram(projCorrStrp(:,c),0.8:0.001:1,'facecolor',vcolors{c},'edgecolor',vcolors{c},'orientation','horizontal','normalization','probability');
end

ylabel('same/other correlation');
xlabel('frequency');
legend('odor','meaning','value');
bootp=(sum(sum(projCorrStrp(:,1)>=projCorrStrp(:,3)'))+1)/(numstraps^2+1);
text(0.05,.97,sprintf('p = %g',round(bootp,3,'significant')));
bootp=(sum(sum(projCorrStrp(:,2)>=projCorrStrp(:,3)'))+1)/(numstraps^2+1);
text(0.15,.97,sprintf('p = %g',round(bootp,3,'significant')));
bootp=(sum(sum(projCorrStrp(:,1)>=projCorrStrp(:,2)'))+1)/(numstraps^2+1);
text(0.2,.93,sprintf('p = %g',round(bootp,3,'significant')));
ylim([0.9 1]);
yticks([0.9 0.95 1]);

%% perform value analysis with trial history
runglm=false;

folds=4;
varExpValueHist=NaN(totalNeurons,1);
%kernels=NaN(totalNeurons,31);
binsPerCue=diff(windows{1})/binSize;

opts.alpha=0.5;
glmopts=glmnetSet(opts);

if runglm
NN=0;
sessInd=find(includedSessions);
tic
for NS=1:sum(includedSessions)
    A = Xall3{NS}(:,[1:binsPerCue size(Xall3{NS},2)-5:size(Xall3{NS},2)]);
    trains = getfolds(A,folds,binsPerTrial);
    NN=NNstart;
    for neuron=1:size(Yall{NS},2)
        NN=NN+1;
        y=Yall{NS}(:,neuron);
        %cross-validated variance explained
        pred=NaN(size(y,1),1);
        for fold=1:folds
            train=trains{fold};
            test=train==0;
            %kernel=lassoglm(A(train,:),y(train),'normal','alpha',0.5,'lambda',thisLambda);
            fit = glmnet(A(train,:),y(train),[],glmopts);
            prediction = glmnetPredict(fit,A(test,:));
            pred(test)=prediction(:,end);
        end
        varExpValueHist(NN,1) = 1- var(y-pred)/var(y);
        kernel = fit.beta(:,end);
        kernels(NN,:)=kernel;
    end    
    NNstart=NN;
    fprintf('Session #%d \n',NS);
end
toc
save('valueGLMFile.mat','varExp','varExpValueSh','varExpValueNew','varExpValueHist','kernels','lambdaVal','predF','-v7.3');

end


%% perform regression for shuffling trial values for value cells
alreadyRun=true;

if ~alreadyRun
tic

shuffles=1000; %how many times to shuffle trial values within cue
folds=4;
NS=0;
varExpValueTrlSh=NaN(totalNeurons,shuffles);
binsPerCue=diff(windows{1})/binSize;
NN=0;

for session=1:sum(includedSessions)
    
    NS=NS+1;
    
    A = Xall{NS};
    trains = getfolds(A,folds,binsPerTrial);

    %get original cue order for shuffle
    firstBins={};
    cueOrder=[1 3 5 2 4 6];
    for cue=1:6
        firstBins{cue}=find(A(:,1+(cue-1)*binsPerCue));
        firstBins{cue}=[firstBins{cue} ones(length(firstBins{cue}),1)*cueOrder(cue)];
    end
    alltrials=cat(1,firstBins{:});
    [sortedtrials,ordered]=sort(alltrials(:,1));
    originalOrder=alltrials(ordered,2);
    
    %make new shuffled predictor matrices
    A3Sh={};
    for perm=1:size(varExpValueTrlSh,2)
        
        newOrder=originalOrder;
        newValue=NaN(length(alltrials),1);
        for cue=1:6
            originalValues=trialValues{session}(originalOrder==cue);
            newValue(newOrder==cue)=originalValues(randsample(sum(originalOrder==cue),sum(newOrder==cue),'false'));
        end
        
        A3Sh{perm}=Xall3{session}(:,[1:binsPerCue size(Xall3{NS},2)-5:size(Xall3{NS},2)]);
        for bin=1:binsPerCue
            A3Sh{perm}(sortedtrials(:,1)-1+bin,bin)=newValue;
        end
    end
    
    for neuron=1:size(Yall{NS},2)
        NN=NN+1;
        if ismember(category(NN,1),[1 2])
            
            y=Yall{NS}(:,neuron);
            
            parfor perm=1:size(varExpValueTrlSh,2)
                A=A3Sh{perm};
                %cross-validated variance explained
                pred=NaN(size(y,1),1);
                for fold=1:folds
                    train=trains{fold};
                    test=train==0;
                    %kernel=lassoglm(A(train,:),y(train),'normal','alpha',0.5,'lambda',thisLambda);
                    fit = glmnet(A(train,:),y(train),[],glmopts);
                    prediction = glmnetPredict(fit,A(test,:));
                    pred(test)=prediction(:,end);
                end
                varExpValueTrlSh(NN,perm) = 1- var(y-pred)/var(y);
                
                
                
            end
            
        end
    end
    fprintf('Session #%d \n',NS);
    
end
save('valueGLMFile.mat','varExp','varExpValueSh','varExpValueNew','varExpValueHist','kernels','varExpValueTrlSh','predF','-v7.3');

toc
end
%% history coding figure
figure;
%categories for subselection

sel = ismember(category,[1:8]); %which neurons are allowed
%usedModels=varExpValueNew(:,91:end);
models2compare = [varExpValueHist varExpValueNew(:,91:end)];
%models2compare = [usedModels(:,valRank(1:17))];
[~,bestModelH]=max(models2compare,[],2);
overShuffle=sum(varExpValueHist>varExpValueTrlSh,2)/size(varExpValueTrlSh,2);
history=bestModelH==1 & overShuffle>0.95 & sel;
neuronRegionOlf=neuronRegionAdj;
neuronRegionOlf(ismember(neuronRegionAdj,{'EPd','PIR'}))={'OLF'};
allReg=unique(neuronRegionOlf);

regions={'ALM','ACA','FRP','PL','ILA','ORB','DP','TTd','AON'};
reggrps=[2:4];
neuronGroupAdj = neuronGroup + 1;
neuronGroupAdj(neuronGroupAdj==5)=1;

% regions={'ACA','FRP','PL','ILA','ORB','CP','ACB'};
% reggrps=[2 4];

allReg=unique(neuronRegionOlf);
regInd=NaN(length(regions),1);
order=[9 10 5 2 3 1 7 12 8 11 6 13 4];
allReg=allReg(order);
for reg=1:length(regions)
    regInd(reg,1)=find(ismember(allReg,regions(reg)));
end

pVals=NaN(36,1);

for ct=1
    %region
    barData=zeros(2,5);
    subplot(2,2,1);
    hold on;
    sel=history;
    %othervalue=ismember(category,1)&~sel;
    othersel=ismember(category,1:8)&~sel;
    regionnumber=zeros(totalNeurons,1);
    for reg=1:length(allReg)
        regsel=ismember(neuronRegionOlf,allReg(reg));
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel);
        %barData(2,reg)=sum(othervalue&regsel)/sum(regsel); 
        barData(2,reg)=sum(othersel&regsel)/sum(regsel);          
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial','dummyvarcoding','reference');
    [ypred,ypredCI]=predict(lm,'conditional',false);
    
    Hall=[];
%     Htemp = [1 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0];
%     for reg1 = 1:length(regions)
%         H=Htemp;
%         H(2,reg1+4)=1;
%         Hall=cat(1,Hall,H);        
%     end   
    Htemp = [0 0 0 0 0 0 0 0 0 0 0 0 0];
    for reg1 = 1:length(regions)-1
        for reg2 = reg1+1:length(regions)
            H=Htemp;
            H(reg1+4)=1;
            H(reg2+4)=-1;
            Hall=cat(1,Hall,H);
        end
    end
    greater=zeros(length(regions));
    for reg1 = 1:length(regions)
        for reg2 = 1:length(regions)
            H=Htemp;
            H(reg1+4)=1;
            H(reg2+4)=-1;
            greater(reg1,reg2)=coefTest(lm,H);
        end
    end    
    
    %Hall=unique(Hall,'rows');
    
    for p=1:length(pVals)
        pVals(p,ct) = coefTest(lm,Hall(p,:));
    end    
    
    upperCI=[];
    lowerCI=[];
    estmean=[];
    for reg=1:length(regions)
        regsel=find(ismember(neuronRegionOlf,regions(reg)),1);
        estmean(1,reg)=ypred(regsel);
        upperCI(1,reg)=ypredCI(regsel,2);
        lowerCI(1,reg)=ypredCI(regsel,1);
    end    
 
    barData=barData(:,regInd);
    b=bar(barData','stacked');
    b(2).FaceColor=[0.85 0.85 0.85]; %meaning
    b(1).FaceColor=[0.4 0.9 1]; %value
   
    
    errorbar(1:length(regions),estmean(1,:),estmean(1,:)-lowerCI(1,:),upperCI(1,:)-estmean(1,:),'o','linewidth',0.75,'color','k');
    
    xlim([0.5 length(regions)+0.5]);
    %legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
    xticks(1:length(regions));
    xticklabels(regions);
    xtickangle(45);
    ylim([0 0.55]);
    
    
    %make heatmap
    target=[0.4 0.9 1];
    customMap=NaN(100,3);
    for i=1:length(customMap)
        customMap(i,:)=target*i/(length(customMap));
    end
    customMap=cat(1,[1 1 1],customMap);
    
    %by region
    %greater=zeros(length(regions));
    frac=zeros(length(regions));
    for r1=1:length(regions)
        for r2=1:length(regions)
            frac(r1,r2)=barData(1,r1)-barData(1,r2);
            %greater(r1,r2)=lowerCI(1,r1)>upperCI(1,r2);
        end        
    end
    s=subplot(2,2,2);
    pcutoff=0.05/36;
    imagesc(frac .* (greater<pcutoff & frac>0), [0 ceil(max(frac,[],'all')*20)/20]);
    
    colormap(s,customMap);
    cb=colorbar;
    set(cb,'ticks',[0:0.05:0.3]);
    xticks(1:length(regions));
    xticklabels(regions);
    xtickangle(45);
    yticks(1:length(regions));
    yticklabels(regions);
    %title(titles{ct});
    if ct==1
        ylabel('increase in this region...');
        xlabel('over this region');
    end    
    
    %region group
    barData=zeros(2,3);
    subplot(2,2,3);
    hold on;
    regionnumber=zeros(totalNeurons,1);
    for reg=1:4
        regsel=ismember(neuronGroupAdj,reg);
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel);
        %barData(2,reg)=sum(othervalue&regsel)/sum(regsel); 
        barData(2,reg)=sum(othersel&regsel)/sum(regsel);              
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial');
    [ypred,ypredCI]=predict(lm,'conditional',false);
    
    barData=barData(:,reggrps);    
    b=bar(barData','stacked');
    %b(2).FaceColor=[0 0.7 1];
    b(2).FaceColor=[0.85 0.85 0.85]; 
    b(1).FaceColor=[0.4 0.9 1];
    
    upperCI=[];
    lowerCI=[];
    estmean=[];
    for reg=1:length(reggrps)
        regsel=find(ismember(neuronGroupAdj,reggrps(reg)),1);
        estmean(1,reg)=ypred(regsel);
        upperCI(1,reg)=ypredCI(regsel,2);
        lowerCI(1,reg)=ypredCI(regsel,1);
    end
    
    errorbar(1:length(reggrps),estmean(1,:),estmean(1,:)-lowerCI(1,:),upperCI(1,:)-estmean(1,:),'o','linewidth',0.75,'color','k');
    
    xlim([0.5 7+0.5]);
    %legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
    xticks(1:length(reggrps));
    xticklabels(regionGroupNames(reggrps-1));
    xtickangle(45);
    ylim([0 0.55]);

    
     Htemp = [0 0 0 0];
    greater=zeros(length(reggrps));
    for reg1 = 1:length(reggrps)
        for reg2 = 1:length(reggrps)
            H=Htemp;
            H(reg1+1)=1;
            H(reg2+1)=-1;
            greater(reg1,reg2)=coefTest(lm,H);
        end
    end     
    
    %by region
    %greater=zeros(length(reggrps));
    frac=zeros(length(reggrps));
    for r1=1:length(reggrps)
        for r2=1:length(reggrps)
            frac(r1,r2)=barData(1,r1)-barData(1,r2);
            %greater(r1,r2)=lowerCI(1,r1)>upperCI(1,r2);
        end        
    end
    s=subplot(2,2,4);
    pcutoff=0.05/3;
    imagesc(frac .* (greater<pcutoff & frac>0), [0 ceil(max(frac,[],'all')*20)/20]);
    
    colormap(s,customMap);
    cb=colorbar;
    set(cb,'ticks',[0:0.05:0.3]);
    xticks(1:length(reggrps));
    xticklabels(regionGroupNames(reggrps-1));
    xtickangle(45);
    yticks(1:length(reggrps));
    yticklabels(regionGroupNames(reggrps-1));
    %title(titles{ct});
    if ct==1
        ylabel('increase in this region...');
        xlabel('over this region');
    end        
    
end

%% history value models schematic
examplePerms=[1 1 1];
exampleNeurons=[31 37 63];
%31 37 63
plotWindow=[0 2.5];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);

binTimesGLMTrl = analysisWindow(1)+binSize/2:binSize:analysisWindow(2)-binSize/2;
plotBinsGLM=binTimesGLMTrl>=plotWindow(1) & binTimesGLMTrl<=plotWindow(2);

cueValues=[1 1 0.5 0.5 0 0];
valcolors{1,1}=[0.2 0.2 0.2];
valcolors{2,1}=[0.5 0.5 0.5];
valcolors{3,1}=[0.35 0.05 0.35];
valcolors{4,1}=[0.4 0.1 0.4];
valcolors{5,1}=[0.5 0.15 0.5];
valcolors{6,1}=[0.6 0.2 0.6];
valcolors{7,1}=[0.7 0.25 0.7];
valcolors{8,1}=[0.8 0.3 0.8];
valcolors{9,1}=[0.1 0.6 0.2];
valcolors{10,1}=[0.2 0.8 0.3];

folds=4;

neuronNumber=[1:length(neuronSession)]';
%get neuron number from order in value category
entries=find(history);
for n=1:length(exampleNeurons)
    nn=entries(exampleNeurons(n));
    ta=[];
    subplot(length(exampleNeurons),3,(n-1)*3+1);
    %subplot(6,6,n);
    hold on;
        for cue=1:10
            activity=PSTHv{cue,1}(nn,plotBins);
            ta=[ta activity];
            plot(binTimes(plotBins),activity,'linewidth',1,'color',valcolors{cue});
            
        end
        
    plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
    axis([plotWindow floor(min(ta)) ceil(max(ta))])
    if n>1
        xticks([]);
    end
    if n==1
        xticks([0 2.5]);
        xlabel('seconds from odor onset');
        ylabel('z-score');
    end
    title(num2str(exampleNeurons(n)));

    %get model fits
    session=neuronSession(nn);
    sessionNeuron=find(neuronNumber(neuronSession==session)==nn);
    sessionFolder=fullfile(direc,sessionSubject{session},sessionDate{session});
    cueTimes = readNPY(fullfile(sessionFolder,'cues.times.npy'));
    rewprob = readNPY(fullfile(sessionFolder,'cues.rewardProbability.npy'));
    odorset = readNPY(fullfile(sessionFolder,'cues.odorSet.npy')); 
    cueOrder=zeros(length(cueTimes),1);
    cueOrder(rewprob==1&odorset==1)=1;
    cueOrder(rewprob==1&odorset==2)=2;
    cueOrder(rewprob==0.5&odorset==1)=3;
    cueOrder(rewprob==0.5&odorset==2)=4;
    cueOrder(rewprob==0&odorset==1)=5;
    cueOrder(rewprob==0&odorset==2)=6;
    csp1 = cueTimes(rewprob==1&odorset==1);
    csp2 = cueTimes(rewprob==1&odorset==2);
    csf1 = cueTimes(rewprob==0.5&odorset==1);
    csf2 = cueTimes(rewprob==0.5&odorset==2);
    csm1 = cueTimes(rewprob==0&odorset==1);
    csm2 = cueTimes(rewprob==0&odorset==2);    

    A = Xall3{session}(:,[1:binsPerCue size(Xall3{session},2)-5:size(Xall3{session},2)]);
    trains = getfolds(A,folds,binsPerTrial);
 
    for permn=1:length(examplePerms)
        
        if permn==2
            for c=1:length(cueTimes)
                cueValue=cueValues(cueOrder(c));
                dm=cueValue*eye(binsPerCue);
                A(c*binsPerTrial-binsPerCue+1:c*binsPerTrial,1:binsPerCue)=dm;
            end
        end

        if permn==3
            A = Xall3{session}(:,[binsPerCue+1:size(Xall3{session},2)]);
        end        
        
        
        y=Yall{session}(:,sessionNeuron);
        pred=NaN(size(y,1),1);
        for fold=1:folds
            train=trains{fold};
            test=train==0;
            fit = glmnet(A(train,:),y(train),[],glmopts);
            prediction = glmnetPredict(fit,A(test,:));
            pred(test)=prediction(:,end);
        end
        
        trialPrediction=reshape(pred,[binsPerTrial length(pred)/binsPerTrial])';
        
        %value
        cue1 = ismember(round(cueTimes),round([csp1;csp2]));
        cue2 = ismember(round(cueTimes),round([csf1;csf2]));
        selections{1,1}=trialValues{session}<0.05; %CS- low
        selections{2,1}=trialValues{session}>=0.05 & trialValues{session}<0.1; %CS- high
        selections{3,1}=trialValues{session}>=0.2 & trialValues{session}<0.25; %CS50 low
        selections{4,1}=trialValues{session}>=0.25 & trialValues{session}<0.3;
        selections{5,1}=trialValues{session}>=0.3 & trialValues{session}<0.35;
        selections{6,1}=trialValues{session}>=0.35 & trialValues{session}<0.4;
        selections{7,1}=trialValues{session}>=0.4 & trialValues{session}<0.45;
        selections{8,1}=trialValues{session}>=0.45 & cue2; %CS50 high
        selections{9,1}=trialValues{session}<0.5 & cue1; %CS+ low
        selections{10,1}=trialValues{session}>=0.5; %CS+ high
        for condition=1:length(selections)
            sel = selections{condition};
            examplePSTH{condition,1}=nanmean(trialPrediction(sel,:));
        end
        
        subplot(length(exampleNeurons),4,(n-1)*4+1+permn);
        hold on;
        for cue=1:10
            activity=examplePSTH{cue,1}(plotBinsGLM);
            plot(binTimesGLMTrl(plotBinsGLM),activity,'linewidth',1,'color',valcolors{cue});
        end
        
        plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
        axis([plotWindow floor(min(ta)) ceil(max(ta))])
        if n>1
            xticks([]);
        end
        if n==1
            xticks([0 2.5]);
            title(sprintf('shuffle %d',examplePerms(permn)));
        end
 
        
    end
end  

%% cue dimension value coding


sels={};
sels{1}=history;
sels{2}=category==1 & history==0;
sels{3}=ismember(category,2) & history==0;
sels{4}=ismember(category,3) & history==0;
sels{5}=behavior&~cueK;
cattitles={'value (history)','value','value-like','other cue','lick, cue+lick'};

valcolors{1,1}=[0.2 0.2 0.2];
valcolors{2,1}=[0.5 0.5 0.5];
valcolors{3,1}=[0.35 0.05 0.35];
valcolors{4,1}=[0.4 0.1 0.4];
valcolors{5,1}=[0.5 0.15 0.5];
valcolors{6,1}=[0.6 0.2 0.6];
valcolors{7,1}=[0.7 0.25 0.7];
valcolors{8,1}=[0.8 0.3 0.8];
valcolors{9,1}=[0.1 0.6 0.2];
valcolors{10,1}=[0.2 0.8 0.3];

numstraps=10000;
amin=-0.3;
amax=1;
cueWindow=[0 2.5];
cueBins=binTimes>=cueWindow(1) & binTimes<=cueWindow(2);

meanWindow=[1 2.5];
meanBins=binTimes>=meanWindow(1) & binTimes<=meanWindow(2);


plotWindow=[-3 10];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);
plotBinTimes=binTimes(plotBins);
baseBins=plotBinTimes<0;

%normalize activity from 0-1
dffPSTH6n={};
dffPSTH3n={};
for session=1
    dffPSTH3all=[];
    for condition=1:3
        dffPSTH3all=cat(2,dffPSTH3all,PSTH3{condition,session});
    end
end
for session=1:2
    for condition=1:3
       dffPSTH3n{condition,session}=PSTH3{condition,session}./max(abs(dffPSTH3all),[],2); 
    end
end


%get activity for dimensions
csmActivity=dffPSTH3n{3,1}(:,cueBins);
cspActivity=dffPSTH3n{1,1}(:,cueBins);
csfActivity=dffPSTH3n{2,1}(:,cueBins);

bigBinSize=5;
bigBins=floor(sum(cueBins)/bigBinSize);
bsf=0;
csmActivityB=[];
cspActivityB=[];
csfActivityB=[];
for bb=1:bigBins
    csmActivityB(:,bb)=mean(csmActivity(:,bsf+1:bsf+bigBinSize),2);
    cspActivityB(:,bb)=mean(cspActivity(:,bsf+1:bsf+bigBinSize),2);
    csfActivityB(:,bb)=mean(csfActivity(:,bsf+1:bsf+bigBinSize),2);
    
    bsf=bsf+bigBinSize;
end

[~,PMind]=max(abs(cspActivityB-csmActivityB),[],2,'linear');
[~,PFind]=max(abs(cspActivityB-csfActivityB),[],2,'linear');
[~,FMind]=max(abs(csfActivityB-csmActivityB),[],2,'linear');

differenceVectorPM=cspActivityB(PMind)-csmActivityB(PMind);
differenceVectorPF=cspActivityB(PFind)-csfActivityB(PFind);
differenceVectorFM=csfActivityB(FMind)-csmActivityB(FMind);


sproj={};
oproj={};
disttime={};



%for getting 0 and 1, use z-score
csmActivity=PSTH3{3,1}(:,cueBins);
cspActivity=PSTH3{1,1}(:,cueBins);
csfActivity=PSTH3{2,1}(:,cueBins);

bsf=0;
for bb=1:bigBins
    csmActivityB(:,bb)=mean(csmActivity(:,bsf+1:bsf+bigBinSize),2);
    cspActivityB(:,bb)=mean(cspActivity(:,bsf+1:bsf+bigBinSize),2);
    csfActivityB(:,bb)=mean(csfActivity(:,bsf+1:bsf+bigBinSize),2);
    
    bsf=bsf+bigBinSize;
end

bootpp={};
valSlope=NaN(numstraps,length(sels));
for cl=1:length(sels)
    sel=sels{cl};
    
    PF1=differenceVectorPF(sel)' * cspActivityB(PFind(sel));
    PF0=differenceVectorPF(sel)' * csfActivityB(PFind(sel));
    PM1=differenceVectorPM(sel)' * cspActivityB(PMind(sel));
    PM0=differenceVectorPM(sel)' * csmActivityB(PMind(sel));
    FM1=differenceVectorFM(sel)' * csfActivityB(FMind(sel));
    FM0=differenceVectorFM(sel)' * csmActivityB(FMind(sel));

    c=0;
    for cue=1:10
        c=c+1;
        sproj{c,1}=(differenceVectorPM(sel)' * PSTHv{cue,1}(sel,:) - PM0) ./ (PM1 - PM0);
        sproj{c,2}=(differenceVectorPF(sel)' * PSTHv{cue,1}(sel,:) - PF0) ./ (PF1 - PF0);
        sproj{c,3}=(differenceVectorFM(sel)' * PSTHv{cue,1}(sel,:) - FM0) ./ (FM1 - FM0);
    end

    
    subplot(3,5,5+cl);
    hold on;
    for cue=1:c
        plot(binTimes(plotBins),sproj{cue,1}(plotBins),'linewidth',1,'color',valcolors{cue});

    end
    axis([-0.5 2.5 amin amax]);
    title(cattitles{cl});
    if cl==1
        
            ylabel('CS- ---> CS+');
        xlabel('seconds from odor onset');
    end
    xticks([0 1 2]);
    yticks([0 1]);
        plot([0 0],[-3 12],':','color','k','linewidth',0.5);
    plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
    
    for cue=1:6
        positionFull(cl,cue)=mean(sproj{cue,1}(meanBins));
    end
     
     %bootstrap a statistic
     incNeur=find(sel);
     position=NaN(numstraps,9);

     bootpp{cl}=NaN(length(sels));
     for strap=1:numstraps
         ns=incNeur(randsample(length(incNeur),length(incNeur),'true'));
         PF1=differenceVectorPF(ns)' * cspActivityB(PFind(ns));
         PF0=differenceVectorPF(ns)' * csfActivityB(PFind(ns));
         PM1=differenceVectorPM(ns)' * cspActivityB(PMind(ns));
         PM0=differenceVectorPM(ns)' * csmActivityB(PMind(ns));
         FM1=differenceVectorFM(ns)' * csfActivityB(FMind(ns));
         FM0=differenceVectorFM(ns)' * csmActivityB(FMind(ns));
         
         
         
         c=0;
         for cue=1:10
             c=c+1;
             sproj{c,1}=(differenceVectorPM(ns)' * PSTHv{cue,1}(ns,meanBins) - PM0) ./ (PM1 - PM0);
             sproj{c,2}=(differenceVectorPF(ns)' * PSTHv{cue,1}(ns,meanBins) - PM0) ./ (PM1 - PM0);
             sproj{c,3}=(differenceVectorFM(ns)' * PSTHv{cue,1}(ns,meanBins) - FM0) ./ (FM1 - FM0);
         end

         
         
         
         for cue=1:c
             position(strap,cue)=mean(sproj{cue,1});
         end
         
        b=[(((.225:0.05:.475)-0.05)/0.45);ones(1,6)]'\position(strap,3:8)';
        valSlope(strap,cl)=b(1);

     end
     
     bootpp{cl}=NaN(c);
     for con1=1:c
         for con2=1:c
             bootpp{cl}(con1,con2)=(sum(sum(position(:,con1)<=position(:,con2)'))+1)/(numstraps^2+1);
         end
     end
     
     
     subplot(3,5,10+cl);
     hold on;
     for cue=1:c
         errorbar(cue,mean(position(:,cue)),std(position(:,cue)),'o','color',valcolors{cue}); %'markerfacecolor',valcolors{cue},
     end
     if cl==1
         ylabel('CS-/CS+ axis');
         xlabel('Cue value');
     end
     ylim([-0.1 1]);
     yticks(0:1);
     xticks([]);
     plot([0 11],[0 0],':','color','k','linewidth',0.5);

     
    
end

catcolors={[0.4 0.9 1],[0 0.7 1],[0.7 0 1],[0.6 0.6 0.6],[0.9 0.2 0]};
subplot(3,5,1);
hold on
for s=1:length(sels)
    errorbar(s,nanmean(valSlope(:,s)),nanstd(valSlope(:,s)),'o','color',catcolors{s},'linewidth',1);
end
bootsp=NaN(length(sels));
for con1=1:length(sels)
    for con2=1:length(sels)
        bootsp(con1,con2)=(sum(sum(valSlope(:,con1)<=valSlope(:,con2)'))+1)/(numstraps^2+1);
    end
end
xticks(1:s);
yticks(0:0.5:1);
ylim([0 1]);
xlim([0.5 s+0.5]);
ylabel('CS50 value slope');
basep=(sum(valSlope<=0,1)+1)/(numstraps+1);

subplot(3,4,2);
hold on
plot([-0.1 1],[-0.1 1],':','linewidth',0.75,'color','k');
scatter(varExpValueHist(history),varExp(history,2),24,catcolors{1},'filled');alpha(0.5);
%histogram(varExpValueHist(history)-varExp(history,2),'facecolor',catcolors{1});
xlabel('Var. exp. with history model');
ylabel('Var. exp. with lick model');
xlim([-0.1 1]);
ylim([-0.1 1]);
xticks(0:0.5:1);
yticks(0:0.5:1);

%% 1000 value shuffles
runglm=true;
%generate all value shuffles

folds=4;
NNstart=0;
%varExpValueNew=NaN(totalNeurons,length(valueShuffles));
%kernels=NaN(totalNeurons,31);
binsPerCue=diff(windows{1})/binSize;

opts.alpha=0.5;
glmopts=glmnetSet(opts);

if runglm
varExpValueShuffles=NaN(totalNeurons,1000);
sessInd=find(includedSessions);
valOrder=[1 1 0.5 0.5 0 0];
tic
for NS=1:sum(includedSessions)
    session=sessInd(incses);
    sessionFolder=fullfile(direc,sessionSubject{session},sessionDate{session});
    cueTimes = readNPY(fullfile(sessionFolder,'cues.times.npy'));
    rewprob = readNPY(fullfile(sessionFolder,'cues.rewardProbability.npy'));
    odorset = readNPY(fullfile(sessionFolder,'cues.odorSet.npy'));
    cueOrder=zeros(length(cueTimes),1);
    cueOrder(rewprob==1&odorset==1)=1;
    cueOrder(rewprob==1&odorset==2)=2;
    cueOrder(rewprob==0.5&odorset==1)=3;
    cueOrder(rewprob==0.5&odorset==2)=4;
    cueOrder(rewprob==0&odorset==1)=5;
    cueOrder(rewprob==0&odorset==2)=6;
    
    A = Xall{NS}(:,[1:binsPerCue size(Xall{NS},2)-5:size(Xall{NS},2)]);
    trains = getfolds(A,folds,binsPerTrial);
    for sh=1:1000
        cueOrderSh=cueOrder(randperm(length(cueOrder)),1);
        for c=1:length(cueTimes)
            cueValue=valOrder(cueOrderSh(c));
            dm=cueValue*eye(binsPerCue);
            A(c*binsPerTrial-binsPerCue+1:c*binsPerTrial,1:binsPerCue)=dm;
        end
        NN=NNstart;
        for neuron=1:size(Yall{NS},2)
            NN=NN+1;
            if category(NN)==1
            y=Yall{NS}(:,neuron);
            %cross-validated variance explained
            pred=NaN(size(y,1),1);
            for fold=1:folds
                train=trains{fold};
                test=train==0;
                %kernel=lassoglm(A(train,:),y(train),'normal','alpha',0.5,'lambda',thisLambda);
                fit = glmnet(A(train,:),y(train),[],glmopts);
                prediction = glmnetPredict(fit,A(test,:));
                pred(test)=prediction(:,end);
            end
            varExpValueShuffles(NN,sh) = 1- var(y-pred)/var(y);
            end

        end           
    end
    fprintf('Session #%d \n',NS);
    NNstart=NN;
end
toc
save('valueGLMFile.mat','varExp','varExpValueSh','varExpValueNew','varExpValueShuffles','varExpValueHist','kernels','varExpValueTrlSh','predF','-v7.3');

end

%% how many neurons better than chance fitting of value?

overShuffle=sum(varExpValueShuffles<=varExpValueNew(:,1),2)/size(varExpValueShuffles,2);
overShuffle=overShuffle(category==1);
%% get proportions of value and meaning and compare between regions for striatum
regions={'ACA','FRP','PL','ILA','ORB','CP','ACB'};


%fit all together, allowing random effect of session to have more power
catcolors={[0 0.7 1],[0.7 0 1],[0.6 0.6 0.6]};
titles={'value','value-like','untuned'};
cats={1,2,3};
for ct=1:3
    %region
    barData=zeros(2,5);
    subplot(4,3,ct);
    hold on;
    sel=ismember(category,cats{ct});
    othersel=ismember(category,1:8)&~sel;
    regionnumber=zeros(totalNeurons,1);
    for reg=1:length(regions)
        regsel=ismember(neuronRegionOlf,regions(reg));
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel);
        barData(2,reg)=sum(othersel&regsel)/sum(regsel);   
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial');
    [ypred,ypredCI]=predict(lm,'conditional',false);
    
    b=bar(barData','stacked');
    %b(3).FaceColor=[0.4 0 0.5]; %odor
    b(2).FaceColor=[0.85 0.85 0.85]; %meaning
    b(1).FaceColor=catcolors{ct}; %value
    %b(4).FaceColor=[0 0 0.2]; %odor
    % b(5).FaceColor=[0.9 0.2 0];
    %b(3).FaceColor=[1 1 1];
    upperCI=[];
    lowerCI=[];
    estmean=[];
    for reg=1:length(regions)
        regsel=find(ismember(neuronRegionOlf(regionnumber>0),regions(reg)),1);
        estmean(1,reg)=ypred(regsel);
        upperCI(1,reg)=ypredCI(regsel,2);
        lowerCI(1,reg)=ypredCI(regsel,1);
    end
    
    errorbar(1:length(regions),estmean(1,:),estmean(1,:)-lowerCI(1,:),upperCI(1,:)-estmean(1,:),'o','linewidth',0.75,'color','k');
    
    xlim([0.5 length(regions)+0.5]);
    %legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
    xticks(1:length(regions));
    xticklabels(regions);
    xtickangle(45);
    ylim([0 0.55]);
    
    
    %make heatmap
    target=catcolors{ct};
    customMap=NaN(100,3);
    for i=1:length(customMap)
        customMap(i,:)=target*i/(length(customMap));
    end
    customMap=cat(1,[1 1 1],customMap);
    
    %by region
    greater=zeros(length(regions));
    frac=zeros(length(regions));
    for r1=1:length(regions)
        for r2=1:length(regions)
            frac(r1,r2)=barData(1,r1)-barData(1,r2);
            greater(r1,r2)=lowerCI(1,r1)>upperCI(1,r2);
        end        
    end
    s=subplot(4,3,3+ct);
    imagesc(frac .* greater, [0 ceil(max(frac,[],'all')*20)/20]);
    
    colormap(s,customMap);
    cb=colorbar;
    set(cb,'ticks',[0:0.1:0.3]);
    xticks(1:length(regions));
    xticklabels(regions);
    xtickangle(45);
    yticks(1:length(regions));
    yticklabels(regions);
    title(titles{ct});
    if ct==1
        ylabel('increase in this region...');
        xlabel('over this region');
    end    
    
    %region group
    barData=zeros(2,3);
    subplot(4,3,6+ct);
    hold on;
    sel=ismember(category,cats{ct});
    othersel=ismember(category,1:8)&~sel;    
    regionnumber=zeros(totalNeurons,1);
    regs=[2 4];
    for reg=1:2
        regsel=ismember(neuronGroup,regs(reg));
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel);   
        barData(2,reg)=sum(othersel&regsel)/sum(regsel);           
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial');
    [ypred,ypredCI]=predict(lm,'conditional',false);
    
    b=bar(barData','stacked');
    %b(3).FaceColor=[0.4 0 0.5]; %odor
    b(2).FaceColor=[0.85 0.85 0.85];
    b(1).FaceColor=catcolors{ct}; %value
    %b(4).FaceColor=[0 0 0.2]; %odor
    % b(5).FaceColor=[0.9 0.2 0];
    %b(3).FaceColor=[1 1 1];
    upperCI=[];
    lowerCI=[];
    estmean=[];
    for reg=1:2
        regsel=find(ismember(neuronGroup(regionnumber>0),regs(reg)),1);
        estmean(1,reg)=ypred(regsel);
        upperCI(1,reg)=ypredCI(regsel,2);
        lowerCI(1,reg)=ypredCI(regsel,1);
    end
    
    errorbar(1:2,estmean(1,:),estmean(1,:)-lowerCI(1,:),upperCI(1,:)-estmean(1,:),'o','linewidth',0.75,'color','k');
    
    xlim([0.5 7+0.5]);
    %legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
    xticks(1:3);
    xticklabels(regionGroupNames(regs));
    xtickangle(45);
    ylim([0 0.55]);
    
    
    %make heatmap
    target=catcolors{ct};
    customMap=NaN(100,3);
    for i=1:length(customMap)
        customMap(i,:)=target*i/(length(customMap));
    end
    customMap=cat(1,[1 1 1],customMap);
    
    %by region
    greater=zeros(2);
    frac=zeros(2);
    for r1=1:2
        for r2=1:2
            frac(r1,r2)=barData(1,r1)-barData(1,r2);
            greater(r1,r2)=lowerCI(1,r1)>upperCI(1,r2);
        end        
    end
    s=subplot(4,3,9+ct);
    imagesc(frac .* greater, [0 ceil(max(frac,[],'all')*20)/20]);
    
    colormap(s,customMap);
    cb=colorbar;
    set(cb,'ticks',[0:0.1:0.3]);
    xticks(1:2);
    xticklabels(regionGroupNames(regs));
    xtickangle(45);
    yticks(1:2);
    yticklabels(regionGroupNames(regs));
    title(titles{ct});
    if ct==1
        ylabel('increase in this region...');
        xlabel('over this region');
    end        
    
end


%% 
function A = makeA(discreteUnfixed,continuousPredictorValues)

%discrete unfixed variables
binTimes=discreteUnfixed.binTimes;
binSize=median(diff(binTimes));
binStarts=binTimes-binSize/2;
continuedBin=diff(binTimes)<binSize+binSize/10;
extraBins=0;
for ev=1:length(discreteUnfixed.times)
    extraBins=extraBins+diff(discreteUnfixed.windows{ev})/binSize;
end

A=zeros(length(continuousPredictorValues{1}),round(extraBins));
binsSoFar=0;
for ev=1:length(discreteUnfixed.times)
    vals=discreteUnfixed.values{ev};
    for entry=1:length(discreteUnfixed.times{ev})
        numBinsPerKernel=round(diff(discreteUnfixed.windows{ev})/binSize);
        [closestBinDiff,closestBin]=min(abs(binStarts-(discreteUnfixed.times{ev}(entry)+discreteUnfixed.windows{ev}(1))));
        if closestBinDiff<binSize
            for kk = 1:numBinsPerKernel
                
                A(closestBin+kk-1,kk+binsSoFar)=vals(entry);
                if closestBin+kk-1>length(continuedBin) || continuedBin(closestBin+kk-1)==0
                    break
                end            
            end
        end 
    end
    binsSoFar=binsSoFar+numBinsPerKernel;
end

%add continuous predictors
for ev=1:length(continuousPredictorValues)
    A=cat(2,A,continuousPredictorValues{ev});
end

end

%% prepare predictors for GLM after reduced rank
function trains = getfolds(rA,folds,binsPerTrial)

numTrials=size(rA,1)/binsPerTrial;
blockNum=[];
for trial=1:numTrials
    blockNum=cat(1,blockNum,ceil(trial*ones(binsPerTrial,1)/3));
end

for fold=1:folds
    trains{fold} = mod(blockNum,folds)~=fold-1;
end

end

%%
function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%   Adam Auton, 9th October 2009
if nargin < 1, m = size(get(gcf,'colormap'),1); end
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b]; 
end 

%% nanste function
%standard error, omitting NaNs
function  ste=nanste(dat,dimension)
ste=nanstd(dat,[],dimension)./sqrt(sum(~isnan(dat),dimension));
end