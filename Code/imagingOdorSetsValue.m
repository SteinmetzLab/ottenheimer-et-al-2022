% Script to analyze olfactory conditioning data across day 3 of each odor set
githubDir = 'D:\GitHub';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'steinmetz-et-al-2019')))
direc = 'D:\GitHub\ottenheimer-et-al-2022\Imaging';
addpath(direc);

%% behavior summary
mice={'pl08','pl10','pl11','pl15','pl16'};
sessions={'o1d3','o2d3'};

figure;

%parameters
ymax=6.5;
stimDuration = 1.5;
rewOnset = 2.5;
smoothing=8; %controls standard deviation of smoothing filter

%parameters
binSize = 0.1;
window = [-3 10];
binEdges = window(1):binSize:window(2);
binTimes = window(1)+binSize/2:binSize:window(2)-binSize/2;

plotrange = [-0.5 6];
%smoothing filter for licking PSTH
licksmoothbins=10; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',smoothing);
lickfilterweights=pdf(halfnormal,0:licksmoothbins);

anticipatoryLicks={};
lickPSTH3={};
lickPSTH4={};

ns=10; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',smoothing);
lickfilterweights=pdf(halfnormal,0:licksmoothbins);

alldat=[];
allpreds=[];
allsubj=[];
trialsback=10;
for session=1:length(sessions)
    sn=sessions{session};
    for mouse=1:length(mice)
        mn=mice{mouse};
        sessionName=append(mn,sn);
        load(fullfile(direc,mn,sn,append(sessionName,'events.mat')));
        trialNo = [1:length(cue)]';
        
        trialTbl=table();
        trialTbl.trialNo=trialNo;
        trialTbl.cue1=ismember(round(cue),round(cue1));
        trialTbl.cue2=ismember(round(cue),round(cue2));
        trialTbl.cue3=ismember(round(cue),round(cue3));
        trialTbl.reward=ismember(round(cue),round(cat(1,cue1,cue2r)));
        
        %get licks/bin
        binnedLicks=NaN(length(cue),length(binEdges)-1);
        antLicks=NaN(length(cue),1);
        conLicks=NaN(length(cue),1);
        preLicks=NaN(length(cue),1);
        for trial = 1:length(cue)
            lickTimesRel = (lick - cue(trial));
            binnedLicks(trial,:)=histcounts(lickTimesRel,binEdges)/binSize;
            antLicks(trial,1) = sum(lickTimesRel>0 & lickTimesRel<2.5);
            conLicks(trial,1) = sum(lickTimesRel>3 & lickTimesRel<4.5);
            preLicks(trial,1) = sum(lickTimesRel>-2.5 & lickTimesRel<0);
        end
        
        trialTbl.antLicks=antLicks;
        trialTbl.preLicks=preLicks;
        trialTbl.conLicks=conLicks;
        trialData{mouse,session}=trialTbl;
        
        
        
        dat=(trialTbl.antLicks)/max(trialTbl.antLicks);
        rew=trialTbl.reward;
        
        prevRew=zeros(length(rew),trialsback);
        prevCueRew=zeros(length(rew),trialsback);
        cue2rew=rew(trialTbl.cue2==1);
        
        cue2prev=ones(sum(trialTbl.cue2==1),trialsback)/2;
        for tb=1:trialsback
            cue2prev(tb+1:end,tb)=cue2rew(1:end-tb);
            prevRew(tb+1:end,tb)=rew(1:end-tb);
        end
        prevCueRew(trialTbl.cue1==1,:)=1;
        prevCueRew(trialTbl.cue2==1,:)=cue2prev;
        
        %how many of each cue experienced
        num1=zeros(length(rew),1);
        num1(trialTbl.cue1)=[1:sum(trialTbl.cue1)]/sum(trialTbl.cue1);
        num2=zeros(length(rew),1);
        num2(trialTbl.cue2)=[1:sum(trialTbl.cue2)]/sum(trialTbl.cue2);
        num3=zeros(length(rew),1);
        num3(trialTbl.cue3)=[1:sum(trialTbl.cue3)]/sum(trialTbl.cue3);
        
        preds=[prevRew prevCueRew trialTbl.cue1 trialTbl.cue2 num1 num2 num3];        
        subj=ones(length(rew),1)*mouse;
        allsubj=[allsubj;subj];
        alldat=[alldat;dat];
        allpreds=[allpreds;preds];
        
        
        
        %smooth lick traces
        smoothedLicks=NaN(length(cue),length(binTimes));
        for trial=1:length(cue)
            for l=1:length(binTimes)
                smoothedLicks(trial,l)=sum(binnedLicks(trial,l-min([l-1 licksmoothbins]):l).*fliplr(lickfilterweights(1:min([l licksmoothbins+1]))))/sum(lickfilterweights(1:min([l licksmoothbins+1])));
            end
        end
        
        %3 conditions
        selections={};
        selections{1,1}=cue1; %cs+
        selections{2,1}=cue2; %cs50
        selections{3,1}=cue3; %cs-
        for condition=1:length(selections)
            sel = ismember(round(cue),round(selections{condition}));
            lickPSTH3{condition,1}(mouse,:)=nanmean(smoothedLicks(sel,:));
            anticipatoryLicks{condition,1}(mouse,session)=nanmean(antLicks(sel)) - nanmean(preLicks(sel));
            
        end
        
        %4 conditions
        selections{1,1}=cue1; %cs+
        selections{2,1}=cue2r; %reward+
        selections{3,1}=cue2u; %reward-
        selections{4,1}=cue3; %cs-
        for condition=1:length(selections)
            sel = ismember(round(cue),round(selections{condition}));
            lickPSTH4{condition,1}(mouse,:)=nanmean(smoothedLicks(sel,:));
        end
        
        %lick PSTH based on trial value from coefficient weights
        if exist('trialValues')
            %previous rewards
            tv=trialValues{mouse,session};
            selections{1,1}=tv<0.1;
            selections{2,1}=tv>=0.1 & tv<0.26;
            selections{3,1}=tv>=0.26 & tv<0.34;
            selections{4,1}=tv>=0.34 & tv<0.42;
            selections{5,1}=tv>=0.42 & tv<0.5;
            selections{6,1}=tv>=0.5;
            for condition=1:length(selections)
                sel = selections{condition};
                lickPSTHp{condition,1}(mouse+(session-1)*length(mice),:)=nanmean(smoothedLicks(sel,:));
            end
        end
        
        
    end
    
        %plot PSTHs
    subplot(2,3,session)
    hold on
    
    colors{1,1}=[0.1 0.6 0.2];
    colors{2,1}=[0.4 0.1 0.4];
    colors{3,1}=[0.3 0.3 0.3];
    xvals = find(binTimes>=plotrange(1) & binTimes<=2.5);
    for condition=1:3
        
        %get values
        psth=nanmean(lickPSTH3{condition,1});
        sem=nanste(lickPSTH3{condition,1},1); %calculate standard error of the mean
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
        
        %get values
        psth=nanmean(lickPSTH4{condition,1});
        sem=nanste(lickPSTH4{condition,1},1); %calculate standard error of the mean
        up=psth+sem;
        down=psth-sem;
        
        %plotting
        p{condition}=plot(binTimes(xvals),psth(xvals),'Color',colors{condition,1},'linewidth',1);
        patch([binTimes(xvals),binTimes(xvals(end):-1:xvals(1))],[up(xvals),down(xvals(end):-1:xvals(1))],colors{condition,1},'EdgeColor','none');alpha(0.5);
    end
    
    patch([0 0 stimDuration stimDuration],[0 ymax ymax 0],[0.6 0.3 0],'edgecolor','none');alpha(0.3);
    plot([rewOnset rewOnset],[0 ymax],'color',[0 0.6 0.3]);
    axis([plotrange 0 ymax]);
    xlabel('Seconds from odor delivery');
    ylabel('Licks/s');
    
    if session==1 legend([p{:}],'CS+','CS50(r)','CS50(u)','CS-','location','NW'); end
    title(['Odor set ' num2str(session)]);
    
end

if exist('trialValues')
    subplot(2,3,6);
    hold on;
    colors{1,1}=[0 0 0.2];
    colors{2,1}=[0 0.1 0.4];
    colors{3,1}=[0 0.3 0.7];
    colors{4,1}=[0 0.5 0.9];
    colors{5,1}=[0 0.7 1];
    colors{6,1}=[0.1 0.9 1];
    xvals = find(binTimes>=plotrange(1) & binTimes<=2.5);
    p={};
    for condition=1:6
        
        psthsubj=NaN(length(mice),size(lickPSTHp{condition,1},2));
        for subj=1:length(mice)
            psthsubj(subj,:)=nanmean(lickPSTHp{condition,1}([subj subj+length(mice)],:));
        end
        
        %get values
        psth=nanmean(psthsubj);
        sem=nanste(psthsubj,1); %calculate standard error of the mean
        up=psth+sem;
        down=psth-sem;
        
        %plotting
        p{condition}=plot(binTimes(xvals),psth(xvals),'Color',colors{condition,1},'linewidth',1.25);
        patch([binTimes(xvals),binTimes(xvals(end):-1:xvals(1))],[up(xvals),down(xvals(end):-1:xvals(1))],colors{condition,1},'EdgeColor','none');alpha(0.1);
    end
    
    %patch([0 0 stimDuration stimDuration],[0 ymax ymax 0],[0.6 0.3 0],'edgecolor','none');alpha(0.3);
    axis([0 2.5 0 5]);
    xlabel('seconds from odor onset');
    ylabel('Licks/s');
    
    %legend([p{:}],'0','1','2','3','4','5');
    title('Cue value');
end

rweights=[];
cweights=[];
for s=1:length(mice)
sel=allsubj==s;
mdl=fitlm(allpreds(sel,:),alldat(sel));
rweights(s,:)=mdl.Coefficients.Estimate(2:trialsback+1);
cweights(s,:)=mdl.Coefficients.Estimate(trialsback+2:trialsback*2+1);
end

%allmice together
sel=allsubj>0;
mdl=fitlm(allpreds(sel,:),alldat(sel));
subplot(2,3,4)
hold on;
errorbar(1:trialsback,mdl.Coefficients.Estimate(2:trialsback+1),mdl.Coefficients.SE(2:trialsback+1),'linewidth',1.5,'color','k');
plot(1:trialsback,rweights,'linewidth',0.5);
plot([0 11],[0 0],':','color','k');
ylim([-0.12 0.18]);
xlim([0 trialsback+1]);
xticks([1 5 10]);
ylabel('coefficient');
xlabel('trials back');
title('previous trial rewarded');

subplot(2,3,5)
hold on;
errorbar(1:trialsback,mdl.Coefficients.Estimate(trialsback+2:trialsback*2+1),mdl.Coefficients.SE(trialsback+2:trialsback*2+1),'linewidth',1.5,'color','k');
plot(1:trialsback,cweights,'linewidth',0.5);
plot([0 11],[0 0],':','color','k');
ylim([-0.12 0.18]);
xlim([0 trialsback+1]);
xticks([1 5 10]);

ylabel('coefficient');
xlabel('trials back');
title('previous cue rewarded');

if ~exist('trialValues')
    %get values from weights
    trialValues={};
    for session=1:length(sessions)
        sn=sessions{session};
        for mouse=1:length(mice)
            mn=mice{mouse};
            
            
            trialTbl=trialData{mouse,session};
            
            dat=(trialTbl.antLicks)/max(trialTbl.antLicks);
            rew=trialTbl.reward;
            
            prevRew=zeros(length(rew),trialsback);
            prevCueRew=zeros(length(rew),trialsback);
            cue2rew=rew(trialTbl.cue2==1);
            
            cue2prev=ones(sum(trialTbl.cue2==1),trialsback)/2;
            for tb=1:trialsback
                cue2prev(tb+1:end,tb)=cue2rew(1:end-tb);
                prevRew(tb+1:end,tb)=rew(1:end-tb);
            end
            prevCueRew(trialTbl.cue1==1,:)=1;
            prevCueRew(trialTbl.cue2==1,:)=cue2prev;
            
            
            
            preds=[prevRew prevCueRew trialTbl.cue1 trialTbl.cue2];
            trialValues{mouse,session}=preds * mdl.Coefficients.Estimate(2:end-3);            
            
        end
    end
end



subplot(2,6,6)
hold on;
for session=1:length(sessions)
    errorbar([1 2],mean(anticipatoryLicks{1,1}),nanste(anticipatoryLicks{1,1},1),'linewidth',1.5,'color',[0.1 0.6 0.2]);
    errorbar([1 2],mean(anticipatoryLicks{2,1}),nanste(anticipatoryLicks{2,1},1),'linewidth',1.5,'color',[0.4 0.1 0.4]);
    errorbar([1 2],mean(anticipatoryLicks{3,1}),nanste(anticipatoryLicks{3,1},1),'linewidth',1.5,'color',[0.3 0.3 0.3]);
    plot([1 2] + (session-1)*2.2,anticipatoryLicks{1,1},'color',[0.1 0.6 0.2],'linewidth',0.5); 
    plot([1 2] + (session-1)*2.2,anticipatoryLicks{2,1},'color',[0.4 0.1 0.4],'linewidth',0.5); 
    plot([1 2] + (session-1)*2.2,anticipatoryLicks{3,1},'color',[0.3 0.3 0.3],'linewidth',0.5);
%     errorbar(1,nanmean(mapCorr(:,1)),nanste(mapCorr(:,1),1),'color',[0.2 0.4 0.6],'marker','o','linewidth',1.5);
%     errorbar(2,nanmean(mapCorr(:,2)),nanste(mapCorr(:,2),1),'color',[0.2 0.4 0.6],'marker','o','linewidth',1.5);
%     plot([0.5 2.5],[0 0],':','color','k','linewidth',0.75);
    ylabel('\Delta anticipatory licks');
    plot([0 3],[0 0],':','color','k','linewidth',0.75);
    xticks([1 2]);
    xticklabels({'1','2'});
    xlabel('odor set');
    xlim([0.75 2.25]);
end


%% count neurons
totalNeurons=0;
 
for mouse=1:length(mice)
    mn=mice{mouse};
    %get ROI activity
    load(fullfile(direc,mn,'d3c','Fall.mat'));
    
    % remove doubles and ROIs on edge of FOV
    iscell = processROIs(iscell,stat,F);
    
    numROIs=sum(iscell(:,1));
    totalNeurons=totalNeurons+numROIs;
end


%% roi analysis
neuCoeff=0.7; %coefficient for subtracting neuropil, this is default

mice={'pl08','pl10','pl11','pl15','pl16'};
sessions={'o1d3','o2d3'};

%glm parameters
smoothbinsLickGLM=25; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',20);
filterweightsLickGLM=pdf(halfnormal,0:smoothbinsLickGLM);

%approximate same trial smoothing filter as ephys experiments
smoothbinsTrlEphys=50; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',15);
filterweightsTrlEphys=pdf(halfnormal,0:smoothbinsTrlEphys);
binSizeEphys=0.02;
smoothbinsTrl=15; %number of previous bins used to smooth
binSizeImaging=1/15;
filterweightsTrl=interp1(0:binSizeEphys:binSizeEphys*smoothbinsTrlEphys,filterweightsTrlEphys,0:binSizeImaging:binSizeEphys*smoothbinsTrlEphys);

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
windows{10}=[0 binSize*2]; %block constants

%model with just value for cue
predictornames{1}='cueValue';
predictornames{2}='licks';
predictornames{3}='boutStart';
predictornames{4}='lickRate';
%windows for these predictors
windows3{1}=[0 2.5];
windows3{2}=[-0.3 1];
windows3{3}=[-0.3 3];
windows3{4}=[-numShiftsB*binSize (numShiftsF+1)*binSize]; %there are numShifts*2+1 total for this, so fudging it a bit
windows3{5}=[0 binSize*2]; %block constants

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
% submodelsels{4}=cuePreds|conPreds;
% submodelsels{5}=lickPreds|conPreds;
% submodelsels{6}=rewPreds|conPreds;

Xall=cell(length(mice)*length(sessions),1);
Xall3=cell(length(mice)*length(sessions),1);
Yall=cell(length(mice)*length(sessions),1);
YallSh=cell(length(mice)*length(sessions),1);

dffPSTH3={};
dffPSTH6={};
dffPSTHv={};
dffPSTHf={};
NS=0;
select1=[1 0];
select2=[0 1];
tic
for session=1:length(sessions)
    sn=sessions{session}; 
    NN=0;
    for mouse=1:length(mice)
        NS=NS+1;
        mn=mice{mouse};
        
        sessionName=append(mn,sn);
        load(fullfile(direc,mn,sn,append(sessionName,'events.mat')));

        rewardswithlicks=reward(reward<lick(end));
        lickreward = arrayfun(@(x) lick(find(lick>reward(x),1)),1:length(rewardswithlicks))';
        lickreward = unique(lickreward);
        
        %only use lick reward when at least 400ms after reward
        lrdelay=[];
        for l=1:length(lickreward)
            ldiffs=lickreward(l)-reward;
            lrdelay(l,1)=min(ldiffs(ldiffs>0));
        end
        lickreward(lrdelay<0.4)=[];
        
        %first lick in bout
        ibi=0.5;
        boutStart = lick(cat(1,1,1+find(diff(lick)>ibi)));
        
        
        %%% set up GLM %%%
        binEdgesLickGLM=cat(1,frameTimes(1)-median(diff(frameTimes))/2,frameTimes(1:end-1)+diff(frameTimes)/2,frameTimes(end)+median(diff(frameTimes))/2);
        binnedLicksGLM=histcounts(lick,binEdgesLickGLM) ./ diff(binEdgesLickGLM)';
        lickRateGLM=NaN(length(binnedLicksGLM),1);
        for l=1:length(binnedLicksGLM)
            lickRateGLM(l,1)=sum(binnedLicksGLM(l-min([l-1 smoothbinsLickGLM]):l).*fliplr(filterweightsLickGLM(1:min([l smoothbinsLickGLM+1]))))/sum(filterweightsLickGLM(1:min([l smoothbinsLickGLM+1])));
        end
        
        binnedLicksGLM=NaN(length(cue),length(binEdges)-1);
        binTimesGLM=NaN(length(cue),length(binEdges)-1);
        for trial = 1:length(cue)
            frameTimesRel = (frameTimes - cue(trial));
            for bin=1:length(binTimes)
                binnedLicksGLM(trial,bin)=mean(lickRateGLM(frameTimesRel>=binEdges(bin)&frameTimesRel<binEdges(bin+1)));
                binTimesGLM(trial,bin)=binEdges(bin)+cue(trial)+binSize/2;
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
        GLMbins.(sn){mouse,1} = binTimesGLMincluded(:);

        %make predictor matrix (called A here)
        continuousPredictors=[vector1 ones(length(vector1{1}),1)*select1(session) ones(length(vector1{1}),1)*select2(session)];
        discreteUnfixed.times={unique(cue1*select1(session)),unique(cue2*select1(session)),unique(cue3*select1(session)),unique(cue1*select2(session)),unique(cue2*select2(session)),unique(cue3*select2(session)),lick,boutStart};
        discreteUnfixed.windows=windows(1:8);
        discreteUnfixed.binTimes=binTimesGLMincluded(:);
        discreteUnfixed.values={ones(length(unique(cue1*select1(session))),1),ones(length(unique(cue2*select1(session))),1),ones(length(unique(cue3*select1(session))),1),...
            ones(length(unique(cue1*select2(session))),1),ones(length(unique(cue2*select2(session))),1),ones(length(unique(cue3*select2(session))),1),...
            ones(length(lick),1),ones(length(boutStart),1)};
        A=makeA(discreteUnfixed,continuousPredictors);
        Xall{NS,1}=A;
        
        %make predictor matrix (called A here)
        continuousPredictors=[vector1 ones(length(vector1{1}),1)*select1(session) ones(length(vector1{1}),1)*select2(session)];
        discreteUnfixed.times={cue,lick,boutStart};
        discreteUnfixed.windows=windows3(1:3);
        discreteUnfixed.binTimes=binTimesGLMincluded(:);
        discreteUnfixed.values={trialValues{mouse,session},...
            ones(length(lick),1),ones(length(boutStart),1)};
        A=makeA(discreteUnfixed,continuousPredictors);
        Xall3{NS,1}=A;

        trialNo = [1:length(cue)]';     
        
        %get ROI activity
        load(fullfile(direc,mn,'d3c','Fall.mat'));
        
        % remove doubles and ROIs on edge of FOV
        iscell = processROIs(iscell,stat,F);
        numROIs=sum(iscell(:,1));
        indx=find(iscell(:,1));
        
        
        %determine which frames belong to which session
        if session==1
            firstFrame(mouse,1)=1;
            firstFrame(mouse,2)=length(frameTimes)+2;
            lastFrame(mouse,1)=length(frameTimes);
            lastFrame(mouse,2)=length(F);
        end

        
        for roi=1:numROIs
            NN=NN+1;
            roiNum=indx(roi);
  
                
%             correctedF=F(roiNum,firstFrame(mouse,session):lastFrame(mouse,session))-(neuCoeff * Fneu(roiNum,firstFrame(mouse,session):lastFrame(mouse,session)));
%             medianF=medfilt1(double(correctedF),7);
            
            decSpks=spks(roiNum,firstFrame(mouse,session):lastFrame(mouse,session));
            %decSpks=smooth(decSpks,5)'; %to account for error in estimating deconvolved spike frame
            
            %smooth with the same filter used for ephys
            spikeRateSmooth=NaN(1,length(decSpks));
            for l=1:length(spikeRateSmooth)
                spikeRateSmooth(1,l)=sum(decSpks(1,l-min([l-1 smoothbinsTrl]):l).*fliplr(filterweightsTrl(1:min([l smoothbinsTrl+1]))))/sum(filterweightsTrl(1:min([l smoothbinsTrl+1])));
            end
            
            if session==1 roiMouse{NN,1}=mn; end
            %get roi activity/bin
            binnedActivity=NaN(length(cue),length(binEdges)-1);
            numBins=length(binTimes);
            for trial = 1:length(cue)
                frameTimesRel = (frameTimes - cue(trial));
                [framesPerBin,~,binNumber]=histcounts(frameTimesRel,binEdges);
                includedFrames=binNumber>0;
                binnedActivity(trial,:)=accumarray(binNumber(includedFrames),spikeRateSmooth(includedFrames),[numBins 1]) ./ framesPerBin';
            end
            binnedActivity(isnan(binnedActivity))=0;binnedActivity(binnedActivity==inf)=0; %this only happens if bins go beyond end of session, so extremely rare
%              %convert to dF/F
%              binneddFF = (binnedActivity - median(binnedActivity(:,binTimes<0),2)) / median(F(roiNum,:));

            %to avoid silly z-scores
            neuronStd(NN,1)=std(binnedActivity(:,binTimes<0),0,'all');
            if std(binnedActivity(:,binTimes<0),0,'all')>=1
                binneddFF=(binnedActivity-mean(binnedActivity(:,binTimes<0),'all')) / std(binnedActivity(:,binTimes<0),0,'all');
            else
                binneddFF=(binnedActivity-mean(binnedActivity(:,binTimes<0),'all'));
            end
            
            %3 conditions
            selections={};
            selections{1,1}=cue1; %cs+
            selections{2,1}=cue2; %cs50
            selections{3,1}=cue3; %cs-
            for condition=1:length(selections)
                sel = ismember(round(cue),round(selections{condition}));
                dffPSTH3{condition,session}(NN,:)=nanmean(binneddFF(sel,:));

            end         

%             %4 conditions
%             selections{1,1}=cue1; %cs+
%             selections{2,1}=cue2r; %reward+
%             selections{3,1}=cue2u; %reward-
%             selections{4,1}=cue3; %cs- 
%             for condition=1:length(selections)
%                 sel = ismember(round(cue),round(selections{condition}));
%                 dffPSTH4{condition,session}(NN,:)=nanmean(binneddFF(sel,:));
%             end
            
            %3 conditions, odd and even trials
            selections={};
            selections{1,1}=cue1; %cs+
            selections{2,1}=cue1; %cs+
            selections{3,1}=cue2; %cs50
            selections{4,1}=cue2; %cs50
            selections{5,1}=cue3; %cs-
            selections{6,1}=cue3; %cs-
            for condition=1:length(selections)
                sel = ismember(round(cue),round(selections{condition})) & rem(trialNo,2)==rem(condition,2);
                dffPSTH6{condition,session}(NN,:)=nanmean(binneddFF(sel,:));

            end              
            
            %value all trials
            tv=trialValues{mouse,session};
            selections{1,1}=tv<0.1;
            selections{2,1}=tv>=0.1 & tv<0.26;
            selections{3,1}=tv>=0.26 & tv<0.34;
            selections{4,1}=tv>=0.34 & tv<0.42;
            selections{5,1}=tv>=0.42 & tv<0.5;
            selections{6,1}=tv>=0.5;
            for condition=1:length(selections)
                sel = selections{condition};
                dffPSTHv{condition,1}(NN,:)=nanmean(binneddFF(sel,:));
            end

            %value CS50
            selections{1,1}=tv<0.24;
            selections{2,1}=tv>=0.24 & tv<0.28;
            selections{3,1}=tv>=0.28 & tv<0.32;
            selections{4,1}=tv>=0.32 & tv<0.35;
            selections{5,1}=tv>=0.35 & tv<0.38;
            selections{6,1}=tv>=0.38;
            for condition=1:length(selections)
                sel = selections{condition} & ismember(round(cue),round(cue2));
                dffPSTHf{condition,1}(NN,:)=nanmean(binneddFF(sel,:));
            end            
            
            %%% prepare for GLM %%%%
            includedActivity=binneddFF(:,includedBins)';
            shOrder=randperm(size(binneddFF,1));
            includedActivitySh=binneddFF(shOrder,includedBins)';
            activityVector=includedActivity(:);
            activityVectorSh=includedActivitySh(:);
            %normalize from 0 to 1
%             activityVectorSh=activityVectorSh/std(activityVectorSh);
%             activityVector=activityVector/std(activityVector);
            Yall{NS,1}(:,roi)=activityVector;
            YallSh{NS,1}(:,roi)=activityVectorSh;
            
        end
        fprintf('Session #%d \n',NS);
    end
end
toc


%% get reduced rank predictors
%combine across days for each mouse first
for mouse=1:length(mice)
    XallC{mouse}=[Xall{mouse};Xall{mouse+length(mice)}];
    YallC{mouse}=[Yall{mouse};Yall{mouse+length(mice)}];
    XallC3{mouse}=[Xall3{mouse};Xall3{mouse+length(mice)}];
end

%use combined sessions
[~, bR, R2] = CanonCor2all(YallC, XallC);
[~, bR3, R2] = CanonCor2all(YallC, XallC3);
%% perform regression with reduced rank predictors
components=20;
valueComponents=10;
folds=4;
NS=0;
varExp=NaN(length(roiMouse),5);
predF=cell(length(roiMouse),length(submodelsels)+1);
kernels=NaN(length(roiMouse),length(R2));
lambdaVal=NaN(length(roiMouse),1);
lambda=[0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.5];
binsPerCue=diff(windows{1})/binSize;

%permute all possible combinations of odors for meaning shuffle
allodors=1:6;
perms6=NaN(90,6);
first2=nchoosek(allodors,2);
perm=0;
for f2=1:length(first2)
    remaining=allodors(~ismember(allodors,first2(f2,:)));
    next2=nchoosek(remaining,2);
    for n2=1:length(next2)
        last=remaining(~ismember(remaining,next2(n2,:)));
        perm=perm+1;
        perms6(perm,:)=[first2(f2,:) next2(n2,:) last];
    end
end
perms6all=perms6; 
varExpValueSh=NaN(totalNeurons,4+length(perms6all));

NN=0;
for mouse=1:length(mice)
    session1=mouse;
    session2=mouse+length(mice);
    mn=mice{mouse};
    
    A = [Xall{session1};Xall{session2}];
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
    
    %get original cue order for shuffle
    firstBins=[];
    cueOrder=[1 3 5 2 4 6];
    for cue=1:6
        firstBins{cue}=find(A(:,1+(cue-1)*binsPerCue));
        firstBins{cue}=[firstBins{cue} ones(length(firstBins{cue}),1)*cueOrder(cue)];
    end
    alltrials=cat(1,firstBins{:});
    [sortedtrials,ordered]=sort(alltrials(:,1));
    originalOrder=alltrials(ordered,2);
    combinedValues=[trialValues{mouse,1};trialValues{mouse,2}];
    
    %make new shuffled predictor matrices
    A3Sh={};
    for perm=1:length(perms6all)+1
        if perm>size(perms6all,1)
            newValue=ones(length(alltrials),1);
        else
            newOrder=NaN(length(alltrials),1);
            for cue=1:6
                newOrder(originalOrder==cue)=perms6all(perm,cue);
            end
            
            newValue=NaN(length(alltrials),1);
            for cue=1:6
                originalValues=combinedValues(originalOrder==cue);
                newValue(newOrder==cue)=mean(originalValues(randsample(sum(originalOrder==cue),sum(newOrder==cue),'true')));
            end
        end
        
        A3Sh{perm}=XallC3{mouse};
        for bin=1:binsPerCue
            A3Sh{perm}(sortedtrials(:,1)-1+bin,bin)=newValue;
        end
    end
    
    for roi=1:size(Yall{session1},2)
        NN=NN+1;
        y1=Yall{session1}(:,roi);
        y2=Yall{session2}(:,roi);
        y=[y1;y2];
        
        %cross-validated variance explained
        predLam=NaN(size(y,1),length(lambda));
        for fold=1:folds
            train=trains{fold};
            test=train==0;
            %                 fit = glmnet(rA(train,:), y(train), 'gaussian', opts);
            %                 this_a = glmnetCoef(fit, lambda);
            %                 fitK = this_a(2:end);
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
        predF{NN,1}=pred;
        
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
            predF{NN,1+sub}=pred;
            varExp(NN,1+sub) = 1- var(y-pred)/var(y);
            
%             %variance over time
%             for tb=1:binsPerTrial
%                 tbsel=tb:binsPerTrial:length(y);
%                 varExpT.(sn){NN,1+sub}(tb)=1- var(y(tbsel)-pred(tbsel))/var(y(tbsel));
%             end
        end
        
        
        
        %models with alternate cue
        %value instead of cue identitiy
        A3 = XallC3{mouse};
        rA3 = A3 * bR3(:,1:valueComponents);
        pred=NaN(size(y,1),1);
        for fold=1:folds
            train=trains{fold};
            test=train==0;
            fitK=lassoglm(rA3(train,:),y(train),'normal','alpha',0.5,'lambda',thisLambda);
            kernel=bR3(:,1:valueComponents)*fitK;
            pred(test)=A3(test,:)*kernel;
        end
        varExp(NN,3+sub) = 1- var(y-pred)/var(y);
        
        if varExp(NN,1)-varExp(NN,2)>0.02
        
            for perm=1:size(perms6all,1)+1
                
                rA3 = A3Sh{perm} * bR3(:,1:valueComponents);
                trains = getfolds(rA3,folds,binsPerTrial);
                
                %         if ismember(perm,topModels)
                %            display(perms6all(perm,:));display([combinedValues(1:10) newValue(1:10)]);
                %         end
                
                
                pred=NaN(size(y,1),1);
                for fold=1:folds
                    train=trains{fold};
                    test=train==0;
                    fitK=lassoglm(rA3(train,:),y(train),'normal','alpha',0.5,'lambda',thisLambda);
                    kernel=bR3(:,1:valueComponents)*fitK;
                    pred(test)=A3Sh{perm}(test,:)*kernel;
                end
                varExpValueSh(NN,3+perm) = 1- var(y-pred)/var(y);
            end
        
        end
        
    end
    
    fprintf('Mouse #%d \n',mouse);
    
end


varExpValueSh(:,1)=varExp(:,1);
varExpValueSh(:,2)=varExp(:,2);
varExpValueSh(:,3)=varExp(:,5);




%% value coding analysis 90 with untuned
improvCutoff=0.02;
overallCutoff=0.02;
improvement=varExp(:,1)-varExp;

%get proportions
predictive = varExp(:,1)>overallCutoff; %more than 2% of variance explained
behavior = improvement(:,3)>improvCutoff; %num licks or lick rate
cueK = improvement(:,2)>improvCutoff; %if best model is at least cutoff better than no cue kernels

totalCueVariance=varExp(:,1)-varExp(:,2); %best model compared to no cue model
valueImp=varExp(:,1)-varExpValueSh; %how much better the best model is than value, shuffles, and one kernel
%correct it so you can't be worse than no cues at all (should I do this?)
for n=1:length(valueImp)
    valueImp(n,valueImp(n,:)>totalCueVariance(n))=totalCueVariance(n);
end

%value specific variance
%whatever is accounted for by value kernel
cueVariance=totalCueVariance - valueImp(:,[4:end]);
[~,bestModel]=max(cueVariance,[],2);
bestModel(isnan(cueVariance(:,1)))=NaN;

figure;

cueNames={'+','+','50','50','-','-'};
odorValues=[0.5 0.5 0.37 0.37 0.05 0.05];
vnames={};
shuffVals=NaN(90,6);
for n=1:length(perms6all)
    vnames{n,1}=append(cueNames{perms6all(n,:)==1},'/',cueNames{perms6all(n,:)==2},', ',cueNames{perms6all(n,:)==3},'/',cueNames{perms6all(n,:)==4},', ',cueNames{perms6all(n,:)==5},'/',cueNames{perms6all(n,:)==6});
    for odor=1:6
        shuffVals(n,odor)=odorValues(perms6all(n,:)==odor);
    end
end

subplot(3,1,1);
hold on
frequency=histcounts(bestModel(cueK&~behavior&~isnan(bestModel)),0.5:1:91.5,'normalization','probability');
outliers=frequency>0.05;
topModels=find(outliers);
b=bar([1:max(bestModel)-1 max(bestModel)+1],frequency,'facecolor',[0 0.2 0.4],'linewidth',0.01);
b.FaceColor='flat';
for e=1:length(topModels)
b.CData(topModels(e),:)=[0.7 0 1]; %trial type
end
b.CData(1,:)=[0 0.7 1]; %value
b.CData(end,:)=[0.6 0.6 0.6]; %non-specific



ylabel('fraction of cue cells');

chance=1/90 * sum(bestModel(cueK&~behavior&~isnan(bestModel))<91)/sum(cueK&~behavior&~isnan(bestModel));
plot([0 91],[chance chance],':','color','k','linewidth',0.5);
xlim([0 93]);
ylim([0 0.3]);
xticks(topModels(1:end-1));
yticks(0:0.1:0.3);
%xtickangle(45);
xticklabels(vnames(outliers(1:end-1)));
xtickangle(45);

category=9*ones(length(improvement),1);
category(cueK&~behavior)=8; %cue neurons
%categorize neurons by their best model if it's in the top models
for bm=1:length(topModels)
category(cueK&~behavior&bestModel==topModels(bm))=bm; %cue, odor
end


%% activity of cue cell types and correlations
figure;
map=redblue(256);
CL=[0 1];
cueWindow=[0 2.5];
cueBins=binTimes>=cueWindow(1) & binTimes<=cueWindow(2);

cueCorrMatrix=NaN(12,12,totalNeurons);
cueCorrSets=NaN(6,6,totalNeurons);
setCorr=NaN(totalNeurons,3);

for NN=1:totalNeurons
    cueVector=NaN(sum(cueBins),12);
    cue3Vector=NaN(sum(cueBins)*3,6);
    v=0;
    u=0;
    for os=1:2
        for condition=1:6
            v=v+1;
            cueVector(:,v)=dffPSTH6{condition,os}(NN,cueBins);
        end
        
    end
    %cueVector=cueVector-mean(cueVector,2);
    corrmat=corrcoef(cueVector);
    
    cueCorrMatrix(:,:,NN)=corrmat;
end


%heatmaps
plotWindow=[-0.5 2.5];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);

sortWindow=[0 1.5];
sortBins=binTimes>=sortWindow(1) & binTimes<=sortWindow(2);



colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];


for ct=1
 sel=category==ct;
thisActivity=dffPSTH3{1,1}(sel,:);
cueResp=mean(thisActivity(:,sortBins),2);
[~,sortOrder]=sort(cueResp);
pn=0;
activity=[];
for cue=1:3
    for os=1:2
        
        pn=pn+1;
        
        ax(1)=subplot(2,20,20+(ct-1)*6+pn);
        colormap(ax(1),map);
        hold on;
        
        
        activity = dffPSTH3{cue,os}(sel,:);
        activity = activity(sortOrder,:);
        

        
        imagesc(binTimes(plotBins),[1 size(activity,1)],activity(:,plotBins),[-5 5]);%[min(min(cspActivity)) max(max(cspActivity))]);
        ylim([0.5 size(activity,1)+0.5])
        plot([0 0],[0.5 sum(sel)+0.5],'color',colors{cue},'linewidth',0.75);
        set(gca,'ytick',[]);
        %if pn==1 ylabel('selective odor'); end
        if pn==1 xlabel('seconds from cue'); end
        if pn>1 xticks([]); end
        

        
    end
    
    
    
end
end



pn=0;
activity=[];
sel=ismember(category,2:6);
thisActivity=dffPSTH3{1,1}(sel,:);
cueResp=mean(thisActivity(:,sortBins),2);

regionCodeOpp=8-category;
sortcrit=[regionCodeOpp(sel) cueResp];
[~,sortOrder]=sortrows(sortcrit,[1 2]);
for cue=1:3
    for os=1:2
        
        pn=pn+1;
        
        ax(1)=subplot(2,20,30+(ct-1)*6+pn);
        colormap(ax(1),map);
        hold on;
        
        
        activity = dffPSTH3{cue,os}(sel,:);
        activity = activity(sortOrder,:);
        

        
        imagesc(binTimes(plotBins),[1 size(activity,1)],activity(:,plotBins),[-5 5]);%[min(min(cspActivity)) max(max(cspActivity))]);
        ylim([0.5 size(activity,1)+0.5])
        plot([0 0],[0.5 sum(sel)+0.5],'color',colors{cue},'linewidth',0.75);
        set(gca,'ytick',[]);
        %if pn==1 ylabel('selective odor'); end
        if pn==1 xlabel('seconds from cue'); end
        if pn>1 xticks([]); end
        

        
    end
    
    
    
end


titles={vnames{outliers(1:end-1)},'untuned','odor'};

for ct=1:7
ax(2)=subplot(4,7,7+ct);
sel=category==ct;
imagesc(nanmean(cueCorrMatrix(:,:,sel),3),CL);
%title(titles(ct));
colormap(ax(2),'bone');
xticks([]);
yticks([]);
end

plotWindow=[-0.5 2.5];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);
cn=0;
for ct=1:7
sel=category==ct;
if sum(sel)>0
    cn=cn+1;
activity=[];
for cue=1:3
    for os=1:2   
        activity = [activity dffPSTH3{cue,os}(sel,plotBins)];
    end
end
activity=activity./max(abs(activity),[],2);
[coeff,score,~,~,explained]=pca(activity');
specs={'-','--'};
for comp=1
    cond=0;
    if explained(comp)>5
    subplot(4,7,ct)
    
    hold on
    for cue=1:3
        for os=1:2
            cond=cond+1;
            plot(binTimes(plotBins),score((cond-1)*sum(plotBins)+1:cond*sum(plotBins),comp),specs{os},'color',colors{cue},'linewidth',1);
        end
    end
    ylabel(sprintf('#%g (%g%%)',comp,round(explained(comp),1)));
    yticks([]);
    xticks([0 2]);
    if comp==1 title(titles{cn}); end
    
%     plot([0 0],[min(score(:,comp)) max(score(:,comp))],':','color','k','linewidth',0.5);
%     plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
    
    end
    xlim(cueWindow);
end
end
end

%% remove doubles (from merging) and edge of FOV rois
function iscell = processROIs(iscell,stat,F)
ROIpixels={};
numROIs=sum(iscell(:,1));
idx=find(iscell(:,1));
flag=zeros(numROIs,1);
for roi=1:numROIs
    roiNum=idx(roi);
    ROIpixels{roi,1}=sub2ind([512 512],stat{roiNum}.ypix+1,stat{roiNum}.xpix+1); %add 1 because it's 0 indexed
    if sum(stat{roiNum}.ypix==0 | stat{roiNum}.ypix==511) | sum(stat{roiNum}.xpix==0 | stat{roiNum}.ypix==511)
        flag(roi,1)=1; %remove cells that are on the edge of FOV
    end
    
    if sum(F(roiNum,:)==0)
        flag(roi,1)=1; %remove cells that have a bug with no fluorescence (rare)
    end    
    
end
overlap=[];
for roi=1:numROIs
    otherROIs=1:numROIs;
    otherROIs(roi)=[];
    otherPixels=unique(cat(2,ROIpixels{otherROIs}));
    overlap(roi,1)=sum(ismember(ROIpixels{roi,1},otherPixels))/length(ROIpixels{roi,1});
end
counting=[1:numROIs]';
remove=(overlap==1 & counting<=(numROIs-sum(overlap==1)/3)) | flag==1;
iscell(idx(remove),1)=0;
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
%% nanste function
%standard error, omitting NaNs
function  ste=nanste(dat,dimension)
ste=nanstd(dat,[],dimension)./sqrt(sum(~isnan(dat),dimension));
end