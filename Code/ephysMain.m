% Script to analyze electrophysiology data
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

%region categories
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
                        
                        %create the ALM region label according to chen et al 2017
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

%% region spatial plotting, 3D
%not enough neurons in EPd or PIR, so rename OLF, which already includes other non-region olfactory neurons
neuronRegionOlf=neuronRegionAdj;
neuronRegionOlf(ismember(neuronRegionAdj,{'EPd','PIR'}))={'OLF'};
division=neuronRegionOlf;
regions={'ALM','MOs','ACA','FRP','PL','ILA','ORB','CP','ACB','DP','TTd','AON','OLF'};
regcolors{1,1}=[1 0.8 0];
regcolors{2,1}=[0.8 0.8 0.6];
regcolors{3,1}=[0.9 0.5 0.2];
regcolors{4,1}=[0.5 0.5 0.1];
regcolors{5,1}=[0.6 0.2 0.6];
regcolors{6,1}=[0.5 0.4 0.7];
regcolors{7,1}=[1 0.6 1];
regcolors{8,1}=[0.4 0.4 1];
regcolors{9,1}=[0.1 0.1 0.6];
regcolors{10,1}=[0.3 0.7 0.5];
regcolors{11,1}=[0.1 0.8 0.8];
regcolors{12,1}=[0.3 0.9 0];
regcolors{13,1}=[0 0.3 0.1];

regionNeurons=[];
for reg=1:length(regions)
    regionNeurons(reg,1)=sum(ismember(division,regions(reg)));
end

neuronXYZmm=neuronXYZ/1000;

%all
subplot(1,1,1);
hold on;
for reg=1:length(regions)
    sel=ismember(division,regions(reg));
    offset=(rand(sum(sel),1)*10-5)/1000;
    p1=scatter3(abs(neuronXYZmm(sel,1))+offset,neuronXYZmm(sel,2)+offset,neuronXYZmm(sel,3)+offset,24,regcolors{reg,1},'filled');
    p1.MarkerFaceAlpha=0.5;
end
view(20,20);
xlabel('ML (mm)');
ylabel('AP (mm)');
zlabel('DV (mm)');

regleg={};
for reg=1:length(regions)
    regleg{reg}=sprintf('%s (%d)',regions{reg},regionNeurons(reg));
end
legend(regleg);
%% analyze behavior

subjects=unique(neuronSubject);

%parameters, used for both behavior and ephys
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

anticipatoryLicks={};
lickPSTH3={};
lickPSTH4={};
sessInd=find(includedSessions);
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
    preLicks=NaN(length(cueTimes),1);
    for trial = 1:length(cueTimes)
        lickTimesRel = (lickTimes - cueTimes(trial));
        binnedLicks(trial,:)=histcounts(lickTimesRel,binEdges)/binSize;
        antLicks(trial,1) = sum(lickTimesRel>0 & lickTimesRel<2.5);
        preLicks(trial,1) = sum(lickTimesRel>-2.5 & lickTimesRel<0);
    end
       
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
    
    
end

for os=1:2
%plot PSTHs
subplot(4,3,os)
hold on

colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];

xvals = find(binTimes>=plotrange(1) & binTimes<=2.5);
for condition=1:3
    
    psthsubj=NaN(length(subjects),size(lickPSTH3{condition,os},2));
    for subject=1:length(subjects)
        psthsubj(subject,:)=nanmean(lickPSTH3{condition,os}(ismember(sessionSubject(includedSessions),subjects(subject)),:));
    end
    
    %get values
    psth=nanmean(psthsubj,1);
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
        psthsubj(subject,:)=nanmean(lickPSTH4{condition,os}(ismember(sessionSubject(includedSessions),subjects(subject)),:));
    end
    
    %get values
    psth=nanmean(psthsubj,1);
    sem=nanste(psthsubj,1); %calculate standard error of the mean
    up=psth+sem;
    down=psth-sem;
    
    %plotting
    p{condition}=plot(binTimes(xvals),psth(xvals),'Color',colors{condition,1},'linewidth',1);
    patch([binTimes(xvals),binTimes(xvals(end):-1:xvals(1))],[up(xvals),down(xvals(end):-1:xvals(1))],colors{condition,1},'EdgeColor','none');alpha(0.5);
end

patch([0 0 stimDuration stimDuration],[0 ymax ymax 0],[0.6 0.3 0],'edgecolor','none');alpha(0.3);
plot([rewOnset rewOnset],[0 ymax],'color',[0 0.6 0.3]);
axis([plotrange 0 ymax]);
xlabel('seconds from odor onset');
ylabel('Licks/s');

if os==2 legend([p{:}],'CS+','CS50(r)','CS50(u)','CS-','location','NW'); end
title(['Odor set ' num2str(os)]);
end
   

colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];

subplot(4,6,6);
hold on;
mat4anov=[];
for cue=1:3
    antLicksSubj=NaN(length(subjects),2);
    for os=1:2
        for subject=1:length(subjects)
            antLicksSubj(subject,os)=nanmean(anticipatoryLicks{cue}(ismember(sessionSubject(includedSessions),subjects(subject)),os));
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

for os=1:2
    subplot(4,3,3+os)
    hold on
    for subject=1:length(subjects)
        scatter(anticipatoryLicks{2}(ismember(sessionSubject(includedSessions),subjects(subject)),os),anticipatoryLicks{1}(ismember(sessionSubject(includedSessions),subjects(subject)),os),18,'filled');
    end
    
    plot([-2 10],[-2 10],':','color','k','linewidth',1);
    axis([-2 10 -2 10])
    xticks([0 5 10]);
    yticks([0 5 10]);
    text(5,1,sprintf('odor set %d',os));
    if os==1
        title('mean anticipatory licks');
        xlabel('CS50 trials');
        ylabel('CS+ trials');
    end
    
    data=[anticipatoryLicks{2}(:,os);anticipatoryLicks{1}(:,os)];
    gl=[zeros(length(data)/2,1);ones(length(data)/2,1)];
    sl=[sessionSubject(includedSessions);sessionSubject(includedSessions)];
    [p,tbl,stats]=anovan(data(:),{gl(:),sl(:)},'varnames',{'cue','subject'},'model','interaction','display','off');
    text(5,0,sprintf('p = %g',round(p(1),2,'significant')));

    subplot(4,3,6+os)
    hold on
    for subject=1:length(subjects)
        scatter(anticipatoryLicks{3}(ismember(sessionSubject(includedSessions),subjects(subject)),os),anticipatoryLicks{1}(ismember(sessionSubject(includedSessions),subjects(subject)),os),18,'filled');
    end
    
    plot([-2 10],[-2 10],':','color','k','linewidth',1);
    axis([-2 10 -2 10])
    xticks([0 5 10]);
    yticks([0 5 10]);
    text(5,1,sprintf('odor set %d',os));
    if os==1
        title('mean anticipatory licks');
        xlabel('CS- trials');
        ylabel('CS+ trials');
    end
 
    data=[anticipatoryLicks{3}(:,os);anticipatoryLicks{1}(:,os)];
    gl=[zeros(length(data)/2,1);ones(length(data)/2,1)];
    sl=[sessionSubject(includedSessions);sessionSubject(includedSessions)];
    [p,tbl,stats]=anovan(data(:),{gl(:),sl(:)},'varnames',{'cue','subject'},'model','interaction','display','off');
    text(5,0,sprintf('p = %g',round(p(1),2,'significant')));    
    
    subplot(4,3,9+os)
    hold on
    for subject=1:length(subjects)
        scatter(anticipatoryLicks{3}(ismember(sessionSubject(includedSessions),subjects(subject)),os),anticipatoryLicks{2}(ismember(sessionSubject(includedSessions),subjects(subject)),os),18,'filled');
    end
    
    plot([-2 10],[-2 10],':','color','k','linewidth',1);
    axis([-2 10 -2 10])
    xticks([0 5 10]);
    yticks([0 5 10]);
    text(5,1,sprintf('odor set %d',os));
    if os==1
        title('mean anticipatory licks');
        xlabel('CS- trials');
        ylabel('CS50+ trials');
    end    
    
    data=[anticipatoryLicks{3}(:,os);anticipatoryLicks{2}(:,os)];
    gl=[zeros(length(data)/2,1);ones(length(data)/2,1)];
    sl=[sessionSubject(includedSessions);sessionSubject(includedSessions)];
    [p,tbl,stats]=anovan(data(:),{gl(:),sl(:)},'varnames',{'cue','subject'},'model','interaction','display','off');
    text(5,0,sprintf('p = %g',round(p(1),2,'significant')));
    
end
%% plot example session trials
figure;
subplot(5,1,1);
hold on;
for trl=1:length(csp1)
    plot([csp1(trl) csp1(trl)],[0 1],'linewidth',0.75,'color',colors{1});
end

for trl=1:length(csp1)
    plot([csp2(trl) csp2(trl)],[0 1],'linewidth',0.75,'color',colors{1});
end

for trl=1:length(csp1)
    plot([csf1(trl) csf1(trl)],[-0.5 0.5],'linewidth',0.75,'color',colors{2});
end

for trl=1:length(csp1)
    plot([csf2(trl) csf2(trl)],[-0.5 0.5],'linewidth',0.75,'color',colors{2});
end

for trl=1:length(csp1)
    plot([csm1(trl) csm1(trl)],[-1 0],'linewidth',0.75,'color',colors{3});
end

for trl=1:length(csp1)
    plot([csm2(trl) csm2(trl)],[-1 0],'linewidth',0.75,'color',colors{3});
end

ylim([-1.5 1]);
yticks([]);
xlabel('session time (seconds)');
%% do analysis on all neurons
preRun=true;
analysisFile='analysis20220608.mat';

analysisWindow=[-1 6.5];
binsPerTrial=diff(analysisWindow)/binSize;
predictornames{1}='csp1';
predictornames{2}='csf1';
predictornames{3}='csm1';
predictornames{4}='csp2';
predictornames{5}='csf2';
predictornames{6}='csm2';
predictornames{7}='licks';
predictornames{8}='boutStart';
predictornames{9}='lickreward';
predictornames{10}='lickRate';

%for shifting lick vector in time
%each shift is binSize*2
numShiftsB=3; %making licks happen later than reality
numShiftsF=2; %making licks happen earlier than reality

windows{1}=[0 5];
windows{2}=[0 5];
windows{3}=[0 5];
windows{4}=[0 5];
windows{5}=[0 5];
windows{6}=[0 5];
windows{7}=[-0.3 0.3];
windows{8}=[-0.3 2];
windows{9}=[0 analysisWindow(end)-2.5];
windows{10}=[-numShiftsB*binSize (numShiftsF+1)*binSize]; %there are numShifts*2+1 total for this, so fudging it a bit
windows{11}=[0 binSize*6]; %block constants


binPred=[];
binsofar=0;
for win=1:length(windows)
    winbins=diff(windows{win})/binSize;
    binPred(binsofar+1:binsofar+winbins)=win;
    binsofar=binsofar+winbins;
end
lickPreds=ismember(binPred,[7 8 10]);
cuePreds=ismember(binPred,[1 2 3 4 5 6]);
rewPreds=ismember(binPred,[9]);
conPreds=ismember(binPred,11);

%removing each variable
submodelsels={};
submodelsels{1}=cuePreds==0;
submodelsels{2}=lickPreds==0;
submodelsels{3}=rewPreds==0;

tic

if preRun
    load(analysisFile)
else

Xall=cell(sum(includedSessions),1);
Yall=cell(sum(includedSessions),1);
YallSh=cell(sum(includedSessions),1);
dffPSTH3=cell(3,2);
dffPSTH4=cell(4,2);
dffPSTH6=cell(6,2);

NS=0;

%smoothing filter for licks
smoothbinsLick=25; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',20);
filterweightsLick=pdf(halfnormal,0:smoothbinsLick);

%smoothing filter for spike train
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
    discreteUnfixed.times={csp1,csf1,csm1,csp2,csf2,csm2,lickTimes,boutStart,lickreward};
    discreteUnfixed.windows=windows(1:9);
    discreteUnfixed.binTimes=binTimesGLMincluded(:);
    A=makeA(discreteUnfixed,continuousPredictors);
    Xall{NS,1}=A;
    
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
                            dffPSTH3{condition,1}(NN,:)=nanmean(normActivity(sel,:));
                        end
                        selections{1,1}=csp2; %cs+
                        selections{2,1}=csf2; %cs50
                        selections{3,1}=csm2; %cs-
                        for condition=1:length(selections)
                            sel = ismember(round(cueTimes),round(selections{condition}));
                            dffPSTH3{condition,2}(NN,:)=nanmean(normActivity(sel,:));
                        end
                        
                        %4 conditions
                        selections={};
                        selections{1,1}=csp1; %cs+
                        selections{2,1}=csf1(ismember(round(csf1),round(reward-2.5))); %reward+
                        selections{3,1}=csf1(~ismember(round(csf1),round(reward-2.5))); %reward-
                        selections{4,1}=csm1; %cs-
                        for condition=1:length(selections)
                            sel = ismember(round(cueTimes),round(selections{condition}));
                            dffPSTH4{condition,1}(NN,:)=nanmean(normActivity(sel,:));
                        end
                        selections{1,1}=csp2; %cs+
                        selections{2,1}=csf2(ismember(round(csf2),round(reward-2.5))); %reward+
                        selections{3,1}=csf2(~ismember(round(csf2),round(reward-2.5))); %reward-
                        selections{4,1}=csm2; %cs-
                        for condition=1:length(selections)
                            sel = ismember(round(cueTimes),round(selections{condition}));
                            dffPSTH4{condition,2}(NN,:)=nanmean(normActivity(sel,:));
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
                            dffPSTH6{condition,1}(NN,:)=nanmean(normActivity(sel,:));
                        end
                        selections{1,1}=csp2; %cs+
                        selections{2,1}=csp2; %cs+
                        selections{3,1}=csf2; %cs50
                        selections{4,1}=csf2; %cs50
                        selections{5,1}=csm2; %cs-
                        selections{6,1}=csm2; %cs-
                        for condition=1:length(selections)
                            sel = ismember(round(cueTimes),round(selections{condition})) & rem(trialNo',2)==rem(condition,2);
                            dffPSTH6{condition,2}(NN,:)=nanmean(normActivity(sel,:));
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

save('newAnalysis.mat','dffPSTH3','dffPSTH4','dffPSTH6','Xall','Yall','YallSh','-v7.3');

end

toc

%% get reduced rank predictors

[~, bR, R2] = CanonCor2all(Yall, Xall);

%% perform regression with reduced rank predictors
preRun=true;
glmFile='GLM20220608';

tic

if preRun
    load(glmFile);
else

components=20;
folds=4;
NS=0;
varExp=NaN(totalNeurons,7);
predF=cell(totalNeurons,1);
kernels=NaN(totalNeurons,1);
lambdaVal=NaN(totalNeurons,1);
lambda=[0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.5];
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
        
        for neuron=1:size(Yall{NS},2)
            NN=NN+1;
            y=Yall{NS}(:,neuron)-mean(Yall{NS}(:,neuron));
            
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
            
            %submodels to find unique variance for each variable
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

%% plot heatmaps of 3 main trial types
map=redblue(256);
figure;
colormap(map);
sel=true(size(dffPSTH3{1,1},1),1);

%heatmaps
cspActivity=dffPSTH3{1,1}; %get the firing rates of neurons of interest
cueWindow=[0 1.5];
cueBins=binTimes>=cueWindow(1) & binTimes<=cueWindow(2);
cueResp=mean(cspActivity(sel,cueBins),2);

regions={'ALM','MOs','ACA','FRP','PL','ILA','ORB','CP','ACB','DP','TTd','AON','OLF'};
regcolors{1,1}=[1 0.8 0];
regcolors{2,1}=[0.8 0.8 0.6];
regcolors{3,1}=[0.9 0.5 0.2];
regcolors{4,1}=[0.5 0.5 0.1];
regcolors{5,1}=[0.6 0.2 0.6];
regcolors{6,1}=[0.5 0.4 0.7];
regcolors{7,1}=[1 0.6 1];
regcolors{8,1}=[0.4 0.4 1];
regcolors{9,1}=[0.1 0.1 0.6];
regcolors{10,1}=[0.3 0.7 0.5];
regcolors{11,1}=[0.1 0.8 0.8];
regcolors{12,1}=[0.3 0.9 0];
regcolors{13,1}=[0 0.3 0.1];

neuronRegionOlf=neuronRegionAdj;
neuronRegionOlf(ismember(neuronRegionAdj,{'EPd','PIR'}))={'OLF'};
regionCode=NaN(length(neuronRegionOlf),1);
for NN=1:length(neuronRegionOlf)
    regionCode(NN,1)=find(ismember(regions,neuronRegionOlf(NN)));
end
regionCodeOpp=12-regionCode;
sortcrit=[regionCodeOpp cueResp];
[~,sortOrder]=sortrows(sortcrit(sel,:),[1 2]);

colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];

plotWindow=[-0.5 6];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);
titles={'CS+       ','CS50      ','CS-     '};
sets={'set1','set2'};
pn=0;
for cue=1:3
    for os=1:length(sets)
        sn=sets{os};
        
        activity = dffPSTH3{cue,os}(sel,:);
        activity=activity(sortOrder,:);
        pn=pn+1;
        subplot(1,6,pn);
        hold on;
        
        imagesc(binTimes(plotBins),[1 sum(sel)],activity(:,plotBins),[-3 3]);%[min(min(cspActivity)) max(max(cspActivity))]);
        ylim([0.5 sum(sel)+0.5])
        plot([0 0],[0.5 sum(sel)+0.5],'color',colors{cue},'linewidth',0.75);
        plot([2.5 2.5],[0.5 sum(sel)+0.5],':','color','k','linewidth',0.25);
        set(gca,'ytick',[]);
        
        %colorbar;
        title(titles{cue});
        if pn==1
            regcodes=unique(regionCodeOpp);
            nsofar=0;
            for reg=1:length(regcodes)
                plot([-0.5 -0.5],[nsofar nsofar+sum(regionCodeOpp==regcodes(reg))],'color',regcolors{end+1-reg},'linewidth',1);
                nsofar=nsofar+sum(regionCodeOpp==regcodes(reg));
            end
            xlabel('seconds from odor onset');
        end
    end
    
    
    
end


%% plot cue onsets for each region
figure;

plotWindow=[-1 2.5];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);
xvals=binTimes(plotBins);

neuronRegionOlf=neuronRegionAdj;
neuronRegionOlf(ismember(neuronRegionAdj,{'EPd','PIR'}))={'OLF'};
division=neuronRegionOlf;

regions={'ALM','PL','ORB','DP'};
regcolors{1,1}=[1 0.8 0];
regcolors{2,1}=[0.6 0.2 0.6];
regcolors{3,1}=[1 0.6 1];
regcolors{4,1}=[0.3 0.7 0.5];


colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];
directions=[1 -1];
for condition=1:3
    diractivity = mean(cat(3,dffPSTH6{(condition-1)*2+1,1},dffPSTH6{(condition-1)*2+2,2}),3);
    direction=sign(max(diractivity,[],2)-abs(min(diractivity,[],2)));
    activity = mean(cat(3,dffPSTH6{(condition-1)*2+2,1},dffPSTH6{(condition-1)*2+1,2}),3);
    for d=1:2
        subplot(2,3,condition+(d-1)*3);
        hold on;
        for reg=1:length(regions)
            regsel=ismember(division,regions(reg));
            
            
            %get values
            psth=nanmean(activity(regsel&direction==directions(d),plotBins));
            sem=nanste(activity(regsel&direction==directions(d),plotBins),1); %calculate standard error of the mean
            up=psth+sem;
            down=psth-sem;
            
            %plotting
            plot(binTimes(plotBins),psth,'Color',regcolors{reg,1},'linewidth',0.75);
            patch([xvals,xvals(end:-1:1)],[up,down(end:-1:1)],regcolors{reg,1},'EdgeColor','none');alpha(0.2);
            
            
            
            
        end
        if condition==1 & d==2
            xlabel('seconds from odor onset');
            ylabel('z-score');
        else
            xticks([]);
        end
        
        if condition>1 yticks([]); end
        plot([0 0],[-1 2],'color',colors{condition},'linewidth',0.75);
        plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
        if d==1 axis([plotWindow -0.25 1.2]); end
        if d==2 axis([plotWindow -0.6 0.2]); end
    end



end

%% plot results of GLM analysis
improvCutoff=0.02;
overallCutoff=0.02;

%truncate negative variance explained (just keep it at 0 prediction)
varExpTruncated=varExp;
varExpTruncated(varExpTruncated<0)=0;
improvement=varExpTruncated(:,1)-varExpTruncated;

%get proportions
category=8*ones(length(improvement),1);
predictive = varExp(:,1)>overallCutoff; %more than 5% of variance explained
behavior = improvement(:,3)>improvCutoff; %num licks or lick rate
rewardK = improvement(:,4)>improvCutoff; %static or changing cue
cueK = improvement(:,2)>improvCutoff; %static or changing cue

category(predictive & behavior,1)=3;
category(predictive & cueK,1)=1;
category(predictive & rewardK,1)=5;
category(predictive & rewardK & cueK,1)=6;
category(predictive & behavior & cueK,1)=2;
category(predictive & behavior & rewardK,1)=4;
category(predictive & behavior & cueK & rewardK,1)=7;

%region
neuronRegionOlf=neuronRegionAdj;
neuronRegionOlf(ismember(neuronRegionAdj,{'EPd','PIR'}))={'OLF'};
division=neuronRegionOlf;
regions={'ALM','ACA','FRP','PL','ILA','ORB','DP','TTd','AON'};

uniqueVarReg={};
varianceReg=[];
sess=unique(neuronSession);
for reg=1:length(regions)
    regSel=ismember(division,regions(reg));
    for s=1:length(sess)
        sesssel=neuronSession==sess(s);
        if sum(sesssel&regSel)>0
            varianceReg(s,reg)=mean(varExpTruncated(regSel&sesssel,1));
            for variable=1:3
                uniqueVarReg{variable}(s,reg)=mean(varExpTruncated(regSel&sesssel,1)-varExpTruncated(regSel&sesssel,1+variable));
            end
        else
            varianceReg(s,reg)=NaN;
            for variable=1:3
                uniqueVarReg{variable}(s,reg)=NaN;
            end
        end
    
    end
end

uniqueVarGrp={};
mice=unique(neuronSubject);
for grp=1:3
    grpSel=ismember(neuronGroup,grp);
    for m=1:length(mice)
        mousesel=ismember(neuronSubject,mice(m));
        if sum(mousesel&grpSel)>0
            for variable=1:3
                uniqueVarGrp{variable}(m,grp)=mean(varExpTruncated(grpSel&mousesel,1)-varExpTruncated(grpSel&mousesel,1+variable));
            end
        else
            for variable=1:3
                uniqueVarGrp{variable}(m,grp)=NaN;
            end
        end    
    end
end

figure('position',[100 100 1200 900]);
subplot(3,4,1);
hold on;

expLim=[0 0.5];
aloLim=[0 0.08];
errorbar(1:length(regions),nanmean(varianceReg),nanste(varianceReg,1),'o','color',[0 0 0],'linewidth',1);
for r=1:length(varianceReg)
scatter(1:length(regions),varianceReg(r,:),24,[0.2 0.2 0.2],'X');
end
p=anova1(varianceReg,[],'off');
if p<0.05 text(2,0.9*expLim(2),'*','fontsize',20); end
xticks(1:length(regions));
xticklabels(regions);
xtickangle(45);
ylim(expLim);
xlim([0.5 length(regions)+0.5]);
ylabel('mean variance explained');


vnames={'cues','licks','reward'};
vcolors={[0 0.4 0.9],[0.9 0.2 0],[0.1 0.7 0]};
subplot(3,4,2);
for variable=1:3
    hold on;
    errorbar(-0.2+(variable-1)*0.2+[1:length(regions)],nanmean(uniqueVarReg{variable}),nanste(uniqueVarReg{variable},1),'o','color',vcolors{variable},'linewidth',1);
  
end
xlim([0.5 length(regions)+0.5]);
ylim(aloLim);
ylabel('unique variance explained');
xlabel('region');
xticks([1:length(regions)]);
xticklabels(regions);
xtickangle(45);
legend(vnames{:},'location','northwest');

subplot(3,7,5);
for variable=1:3
    hold on;
    errorbar(1:3,nanmean(uniqueVarGrp{variable}),nanste(uniqueVarGrp{variable},1),'o','color',vcolors{variable},'linewidth',1);

    plot(1:3,uniqueVarGrp{variable},'color',vcolors{variable},'linewidth',0.5);

  
end
xlim([0.5 3.5]);
ylim(aloLim);
ylabel('unique variance explained');
xticks([1:3]);
xticklabels(regionGroupNames(1:3));
xtickangle(45);

data=cat(1,uniqueVarGrp{:});
ml=repmat([1;2;3;4;5],3,3);
vl=cat(1,ones(5,3),2*ones(5,3),3*ones(5,3));
gl=cat(2,ones(15,1),2*ones(15,1),3*ones(15,1));
% [p,tbl,stats]=anovan(data(:),{vl(:),gl(:)},'varnames',{'var','regGroup'},'model','interaction');
% c=multcompare(stats,'dimension',[1 2]);%,'estimate','column');

subplot(3,3,4);
barData=zeros(8,2);
for ct=1:8
    sel=ismember(category,ct);
    for reg=1:length(regions)
        regsel=ismember(division,regions(reg));
        barData(ct,reg)=sum(sel&regsel)/sum(regsel);      
    end
end

b=bar(barData','stacked');
b(1).FaceColor=[0 0.4 0.9]; %cue
b(2).FaceColor=[0.6 0 0.6];
b(3).FaceColor=[0.9 0.2 0]; %lick
b(4).FaceColor=[0.7 0.5 0]; 
b(5).FaceColor=[0.1 0.7 0]; %reward
b(6).FaceColor=[0 0.5 0.7];
b(7).FaceColor=[0.5 0.5 0.5]; %all
b(8).FaceColor=[1 1 1];
xlim([0.5 length(regions)+0.5]);
legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
xlabel('region');
xticks(1:length(regions));
xticklabels(regions);
xtickangle(45);

subplot(3,5,8);
barData=zeros(8,3);
for ct=1:8
    sel=ismember(category,ct);
    for reg=1:3
        regsel=ismember(neuronGroup,reg);
        barData(ct,reg)=sum(sel&regsel)/sum(regsel);
    end
   
end

b=bar(barData','stacked');
b(1).FaceColor=[0 0.4 0.9]; %cue
b(2).FaceColor=[0.6 0 0.6];
b(3).FaceColor=[0.9 0.2 0]; %lick
b(4).FaceColor=[0.7 0.5 0]; 
b(5).FaceColor=[0.1 0.7 0]; %reward
b(6).FaceColor=[0 0.5 0.7];
b(7).FaceColor=[0.5 0.5 0.5]; %all
b(8).FaceColor=[1 1 1];
xlim([0.5 5.5]);
xticks(1:3);
xticklabels(regionGroupNames(1:3));
xtickangle(45);
%% region comparisons of cue, lick, and both cells using generalized linear mixed effects model
figure;
regions={'ALM','ACA','FRP','PL','ILA','ORB','DP','TTd','AON'};
reggrps=[1:3];

allReg=unique(neuronRegionOlf);
regInd=NaN(length(regions),1);
for reg=1:length(regions)
    regInd(reg,1)=find(ismember(allReg,regions(reg)));
end

%fit all together, allowing random effect of session to have more power
catcolors={[0 0.4 0.9],[0.9 0.2 0],[0.6 0 0.6]};
titles={'cues','licks','both'};
cats=[1 3 2];
for ct=1:3
    %region
    barData=zeros(2,5);
    subplot(4,3,ct);
    hold on;
    sel=ismember(category,cats(ct));
    regionnumber=zeros(totalNeurons,1);
    for reg=1:length(allReg)
        regsel=ismember(neuronRegionOlf,allReg(reg));
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel);      
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial');
    [ypred,ypredCI]=predict(lm,'conditional',false);
    
    barData=barData(:,regInd);
    b=bar(barData','stacked');

    b(1).FaceColor=catcolors{ct}; %value

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
    ylim([0 0.5]);
    
    
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
    sel=ismember(category,cats(ct));
    regionnumber=zeros(totalNeurons,1);
    for reg=1:4
        regsel=ismember(neuronGroup,reg);
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel);      
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial');
    [ypred,ypredCI]=predict(lm,'conditional',false);
    
    barData=barData(:,reggrps);
    b=bar(barData','stacked');
    b(1).FaceColor=catcolors{ct}; %value
    upperCI=[];
    lowerCI=[];
    estmean=[];
    for reg=1:length(reggrps)
        regsel=find(ismember(neuronGroup,reggrps(reg)),1);
        estmean(1,reg)=ypred(regsel);
        upperCI(1,reg)=ypredCI(regsel,2);
        lowerCI(1,reg)=ypredCI(regsel,1);
    end
    
    errorbar(1:length(reggrps),estmean(1,:),estmean(1,:)-lowerCI(1,:),upperCI(1,:)-estmean(1,:),'o','linewidth',0.75,'color','k');
    
    xlim([0.5 7+0.5]);
    %legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
    xticks(1:length(reggrps));
    xticklabels(regionGroupNames(reggrps));
    xtickangle(45);
    ylim([0 0.5]);
    
    
    %make heatmap
    target=catcolors{ct};
    customMap=NaN(100,3);
    for i=1:length(customMap)
        customMap(i,:)=target*i/(length(customMap));
    end
    customMap=cat(1,[1 1 1],customMap);
    
    %by region
    greater=zeros(3);
    frac=zeros(3);
    for r1=1:length(reggrps)
        for r2=1:length(reggrps)
            frac(r1,r2)=barData(1,r1)-barData(1,r2);
            greater(r1,r2)=lowerCI(1,r1)>upperCI(1,r2);
        end        
    end
    s=subplot(4,3,9+ct);
    imagesc(frac .* greater, [0 ceil(max(frac,[],'all')*20)/20]);
    
    colormap(s,customMap);
    cb=colorbar;
    set(cb,'ticks',[0:0.1:0.3]);
    xticks(1:length(reggrps));
    xticklabels(regionGroupNames(reggrps));
    xtickangle(45);
    yticks(1:length(reggrps));
    yticklabels(regionGroupNames(reggrps));
    title(titles{ct});
    if ct==1
        ylabel('increase in this region...');
        xlabel('over this region');
    end        
    
end

%% visualize GLM fitting

%1380
neuron=4050;
%148
trls=148;
ylims=[-3 5];
plotWindow=analysisWindow;
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);

for tr=1:length(trls)
    figure;
    trial=trls(tr);
NS=neuronSession(neuron);
NN=neuron-find(ismember(neuronSession,neuronSession(neuron)),1)+1;
binsToPlot=(trial-1)*binsPerTrial+1:trial*binsPerTrial;
xvals=analysisWindow(1)+binSize/2:binSize:analysisWindow(end)-binSize/2;

subplot(2,2,1);
hold on;
plot(xvals,Yall{NS}(binsToPlot,NN),'linewidth',1,'color','k');
plot(xvals,predF{neuron,1}(binsToPlot),':','linewidth',1,'color','k');

lickts=xvals(Xall{NS}(binsToPlot,304)==1)-0.05;
lickrs=xvals(Xall{NS}(binsToPlot,330)==1)-0.05;
for l=1:length(lickrs)
    plot([lickrs(l) lickrs(l)],ylims,'color',[0.1 0.7 0],'linewidth',1);
end
for l=1:length(lickts)
    plot([lickts(l) lickts(l)],[-3 -2],'color',[0.9 0.2 0],'linewidth',0.75);
end

xlabel('seconds from odor onset');
xlim(analysisWindow);
ylim(ylims);
yticks(-4:4:8);
ylabel('activity (z-score)');

patch([0 0 stimDuration stimDuration],[ylims fliplr(ylims)],[0 0.4 0.9],'edgecolor','none');alpha(0.3);
legend('activity','prediction');

vcolors={[0 0 0],[0 0.4 0.9],[0.9 0.2 0],[0.1 0.7 0]};


subplot(2,2,2);
hold on;
plot(xvals,predF{neuron,1}(binsToPlot),':','linewidth',1,'color','k');
plot(xvals,predF{neuron,2}(binsToPlot),':','linewidth',1,'color',[0 0.4 0.9]);
plot(xvals,predF{neuron,3}(binsToPlot),':','linewidth',1,'color',[0.9 0.2 0]);
plot(xvals,predF{neuron,4}(binsToPlot),':','linewidth',1,'color',[0.1 0.7 0]);

for l=1:length(lickrs)
    plot([lickrs(l) lickrs(l)],ylims,'color',[0.1 0.7 0],'linewidth',1);
end
for l=1:length(lickts)
    plot([lickts(l) lickts(l)],[-3 -2],'color',[0.9 0.2 0],'linewidth',0.75);
end
patch([0 0 stimDuration stimDuration],[ylims fliplr(ylims)],[0 0.4 0.9],'edgecolor','none');alpha(0.3);

legend('prediction',...
    sprintf('w/o cues (%g%%)',round(varExpTruncated(neuron,2),2)*100),...
    sprintf('w/o licks (%g%%)',round(varExpTruncated(neuron,3),2)*100),...
    sprintf('w/o reward (%g%%)',round(varExpTruncated(neuron,4),2)*100));
xlim(analysisWindow);
ylim(ylims);
yticks([]);
xticks([]);

binTimesGLMTrl = analysisWindow(1)+binSize/2:binSize:analysisWindow(2)-binSize/2;
plotBinsGLM=binTimesGLMTrl>=plotWindow(1) & binTimesGLMTrl<=plotWindow(2);

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

    
    subplot(2,2,3);
    hold on
    examplePSTH=[];
    line={};
    line{1}=plot(binTimesGLMTrl(plotBinsGLM),mean([dffPSTH3{1,1}(neuron,plotBins);dffPSTH3{1,2}(neuron,plotBins)]),'color','k','linewidth',1);
    for model=1:4
        
        trialPrediction=reshape(predF{neuron,model},[binsPerTrial length(predF{neuron,model})/binsPerTrial])';
        
        %CS+ prediction
        selections={};
        selections{1,1}=cat(1,csp1,csp2); %cs+
        sel = ismember(round(cueTimes),round(selections{1}));
        examplePSTH(model,:)=nanmean(trialPrediction(sel,:));  
        
        activity=examplePSTH(model,plotBinsGLM);
        line{1+model}=plot(binTimesGLMTrl(plotBinsGLM),activity,':','linewidth',1,'color',vcolors{model});
        
    end
            

    plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
    axis([plotWindow -2 5])

legend('activity',sprintf('prediction (%g%%)',round(varExp(neuron,1),2)*100),...
    sprintf('w/o cues (%g%%)',round(varExpTruncated(neuron,2),2)*100),...
    sprintf('w/o licks (%g%%)',round(varExpTruncated(neuron,3),2)*100),...
    sprintf('w/o reward (%g%%)',round(varExpTruncated(neuron,4),2)*100));

yticks(-4:4:8);

end
%% GLM spatial plotting, 3D

colors={[0 0.4 0.9];... %cue
        [0.9 0.2 0];... %lick
        [0.6 0 0.6];...
        [0.7 0.7 0.7]}; 

%all
neuronXYZmm=neuronXYZ/1000;
cats=[1 3];
for ct=1:2
    subplot(1,3,1+ct);
    hold on;
    sel=category==cats(ct);
    offset=(rand(sum(sel),1)*10-5)/1000;
    p1=scatter3(abs(neuronXYZmm(sel,1))+offset,neuronXYZmm(sel,2)+offset,neuronXYZmm(sel,3)+offset,24,colors{ct},'filled');
    p1.MarkerFaceAlpha=0.2;
view(20,20);
axis([0 3 1 3 -7 -1]);



    xticks([]);
    yticks([]);
    zticks([]);


    subplot(1,3,1);
    hold on;
    offset=(rand(sum(sel),1)*10-5)/1000;
    p1=scatter3(abs(neuronXYZmm(sel,1))+offset,neuronXYZmm(sel,2)+offset,neuronXYZmm(sel,3)+offset,24,colors{ct},'filled');
    p1.MarkerFaceAlpha=0.2;
view(20,20);
axis([0 3 1 3 -7 -1]);
xlabel('ML (mm)');
ylabel('AP (mm)');
zlabel('DV (mm)');
end


%% plot traces of cue, lick, and both cells


colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];

map=redblue(256);
figure('position',[100 100 500 1200]);
colormap(map);
%heatmaps
cueWindow=[0 1.5];
cueBins=binTimes>=cueWindow(1) & binTimes<=cueWindow(2);


baseWindow=[-3 0];
baseBins=binTimes>=baseWindow(1) & binTimes<=baseWindow(2);

plotWindow=[-0.5 6];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);
xvals=binTimes(plotBins);
titles={'CS+','CS50','CS-'};
sets={'set1','set2'};
en=0;

cats=[1 3 2 4];
catnames={'cue cells','lick cells','both cells','others'};
pn=0;
for type=1:3
    catnum=cats(type) ;
    sel=category==catnum;
    cspActivity=dffPSTH3{1}(sel,:); %get the firing rates of neurons of interest
    cueResp=mean(cspActivity(:,cueBins),2);
    
    [xx,sortOrder]=sort(cueResp);
    
    
    
    activity = mean(cat(3,dffPSTH3{1,1}(sel,:),dffPSTH3{1,2}(sel,:)),3);
    activityS=activity(sortOrder,:);
    pn=pn+1;
    subplot(1,6,pn);
    hold on;
    
    imagesc(binTimes(plotBins),[1 size(cspActivity,1)],activityS(:,plotBins),[-5 5]);%[min(min(cspActivity)) max(max(cspActivity))]);
    ylim([0.5 sum(sel)+0.5])
    plot([0 0],[0.5 sum(sel)+0.5],'color',colors{1},'linewidth',0.75);
    plot([2.5 2.5],[0.5 sum(sel)+0.5],':','color','k','linewidth',0.25);
    
    plot([2.5 2.5],[0.5 length(cspActivity)+0.5],'color','k','linewidth',0.5);
    xlim(plotWindow);

        set(gca,'ytick',[]);

    title(sprintf('%s (%d)',catnames{type},sum(sel)));
    xticks([]);
    
end


catcolors={[0 0.4 0.9];... %cue
    
[0.9 0.2 0];... %lick
[0.6 0 0.6];... %both
[0.7 0.7 0.7]};
directions=[1 -1];
for condition=1:3
    diractivity = mean(cat(3,dffPSTH6{(condition-1)*2+1,1},dffPSTH6{(condition-1)*2+2,2}),3);
    direction=sign(max(diractivity,[],2)-abs(min(diractivity,[],2)));
    activity = mean(cat(3,dffPSTH6{(condition-1)*2+2,1},dffPSTH6{(condition-1)*2+1,2}),3);
    for d=1:2
        subplot(2,6,condition+3+(d-1)*6);
        hold on;
        for reg=1:length(cats)
            regsel=ismember(category,cats(reg));
            if cats(reg)==4 regsel=category>3; end
            
            %get values
            psth=nanmean(activity(regsel&direction==directions(d),plotBins));
            sem=nanste(activity(regsel&direction==directions(d),plotBins),1); %calculate standard error of the mean
            up=psth+sem;
            down=psth-sem;
            
            %plotting
            plot(binTimes(plotBins),psth,'Color',catcolors{reg,1},'linewidth',0.75);
            patch([xvals,xvals(end:-1:1)],[up,down(end:-1:1)],catcolors{reg,1},'EdgeColor','none');alpha(0.2);
            
            
            
            
        end
        if condition==1 & d==2
            xlabel('seconds from odor onset');
            ylabel('z-score');
        else
            xticks([]);
        end
        
        if condition>1 yticks([]); end
        plot([2.5 2.5],[-2 4],'color','k','linewidth',0.75);
        patch([0 0 stimDuration stimDuration],[-2 4 4 -2],colors{condition},'edgecolor','none');alpha(0.3);
        plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
        if d==1 axis([plotWindow -0.25 2.25]);yticks(0:1:2); end
        if d==2 axis([plotWindow -1 0.2]);yticks(-1:0.5:0); end
        %if condition==1 legend(regions); end
    end
    
    
    
end




%% region comparisons of cue, lick, and both cells for striatum and pfc
figure;
regions={'ACA','FRP','PL','ILA','ORB','CP','ACB'};
reggrps=[2 4];

allReg=unique(neuronRegionOlf);
regInd=NaN(length(regions),1);
for reg=1:length(regions)
    regInd(reg,1)=find(ismember(allReg,regions(reg)));
end

%fit all together, allowing random effect of session to have more power
catcolors={[0 0.4 0.9],[0.9 0.2 0],[0.6 0 0.6]};
titles={'cues','licks','both'};
cats=[1 3 2];
for ct=1:3
    %region
    barData=zeros(2,5);
    subplot(4,3,ct);
    hold on;
    sel=ismember(category,cats(ct));
    regionnumber=zeros(totalNeurons,1);
    for reg=1:length(allReg)
        regsel=ismember(neuronRegionOlf,allReg(reg));
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel);      
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial');
    [ypred,ypredCI]=predict(lm,'conditional',false);
    
    barData=barData(:,regInd);
    b=bar(barData','stacked');
    b(1).FaceColor=catcolors{ct}; %value
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
    ylim([0 0.5]);
    
    
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
    sel=ismember(category,cats(ct));
    regionnumber=zeros(totalNeurons,1);
    for reg=1:4
        regsel=ismember(neuronGroup,reg);
        regionnumber(regsel)=reg;
        barData(1,reg)=sum(sel&regsel)/sum(regsel);      
    end
    tbl = table(categorical(neuronSession(regionnumber>0)), categorical(regionnumber(regionnumber>0)), sel(regionnumber>0), 'VariableNames', {'session', 'region','val'});
    lm = fitglme(tbl,'val~1+region+(1|session)','distribution','binomial');
    [ypred,ypredCI]=predict(lm,'conditional',false);
    
    barData=barData(:,reggrps);
    b=bar(barData','stacked');
    b(1).FaceColor=catcolors{ct}; %value
    upperCI=[];
    lowerCI=[];
    estmean=[];
    for reg=1:length(reggrps)
        regsel=find(ismember(neuronGroup,reggrps(reg)),1);
        estmean(1,reg)=ypred(regsel);
        upperCI(1,reg)=ypredCI(regsel,2);
        lowerCI(1,reg)=ypredCI(regsel,1);
    end
    
    errorbar(1:length(reggrps),estmean(1,:),estmean(1,:)-lowerCI(1,:),upperCI(1,:)-estmean(1,:),'o','linewidth',0.75,'color','k');
    
    xlim([0.5 7+0.5]);
    %legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
    xticks(1:length(reggrps));
    xticklabels(regionGroupNames(reggrps));
    xtickangle(45);
    ylim([0 0.5]);
    
    
    %make heatmap
    target=catcolors{ct};
    customMap=NaN(100,3);
    for i=1:length(customMap)
        customMap(i,:)=target*i/(length(customMap));
    end
    customMap=cat(1,[1 1 1],customMap);
    
    %by region
    greater=zeros(length(reggrps));
    frac=zeros(length(reggrps));
    for r1=1:length(reggrps)
        for r2=1:length(reggrps)
            frac(r1,r2)=barData(1,r1)-barData(1,r2);
            greater(r1,r2)=lowerCI(1,r1)>upperCI(1,r2);
        end        
    end
    s=subplot(4,3,9+ct);
    imagesc(frac .* greater, [0 ceil(max(frac,[],'all')*20)/20]);
    
    colormap(s,customMap);
    cb=colorbar;
    set(cb,'ticks',[0:0.1:0.3]);
    xticks(1:length(reggrps));
    xticklabels(regionGroupNames(reggrps));
    xtickangle(45);
    yticks(1:length(reggrps));
    yticklabels(regionGroupNames(reggrps));
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
    for entry=1:length(discreteUnfixed.times{ev})
        numBinsPerKernel=round(diff(discreteUnfixed.windows{ev})/binSize);
        [closestBinDiff,closestBin]=min(abs(binStarts-(discreteUnfixed.times{ev}(entry)+discreteUnfixed.windows{ev}(1))));
        if closestBinDiff<binSize
            for kk = 1:numBinsPerKernel
                
                A(closestBin+kk-1,kk+binsSoFar)=1;
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