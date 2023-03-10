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

        end
        

        

        
    end
    
    fprintf('Mouse #%d \n',mouse);
    
end


%% perform value analysis with a ton of shuffles
runglm=true;
%generate all value shuffles
valueShuffles = [];
valueVersion = [];
valueLevels = [1 1 0.5 0.5 0 0;1 1 1 1 0 0;1 1 1 0 0 0;1 1 0 0 0 0;1 0 0 0 0 0;1 1 1 1 1 0;1 1 1 1 1 1];
for vl=1:size(valueLevels,1)
allCombos = perms(valueLevels(vl,:));
allCombos = unique(allCombos,'rows');
valueShuffles = cat(1,valueShuffles,allCombos);
valueVersion = cat(1,valueVersion,vl*ones(size(allCombos,1),1));
end
valueShuffles=fliplr(valueShuffles);

folds=4;
NNstart=0;
varExpValueNew=NaN(totalNeurons,length(valueShuffles));
%kernels=NaN(totalNeurons,31);
binsPerCue=diff(windows{1})/binSize;

opts.alpha=0.5;
glmopts=glmnetSet(opts);

if runglm
cueOrdering=[1 3 5 2 4 6];
tic
for NS=1:length(mice)
    
    cueStarts={};
    for cue=1:6
        cueStarts{cue,1}(:,1)=find(XallC{NS}(:,(cue-1)*binsPerCue+1));
        cueStarts{cue,1}(:,2)=cueOrdering(cue)*ones(length(cueStarts{cue,1}),1);
    end
    
    sortedStarts=sortrows(cat(1,cueStarts{:}));
    cueOrder=sortedStarts(:,2);
    
    A = XallC{NS}(:,[1:binsPerCue size(XallC{NS},2)-1:size(XallC{NS},2)]);
    trains = getfolds(A,folds,binsPerTrial);
    for s=1:length(valueShuffles)
        for c=1:length(cueOrder)
            cueValue=valueShuffles(s,cueOrder(c));
            dm=cueValue*eye(binsPerCue);
            A(c*binsPerTrial-binsPerCue+1:c*binsPerTrial,1:binsPerCue)=dm;
        end
        NN=NNstart;
        for neuron=1:size(Yall{NS},2)
            NN=NN+1;
            y=YallC{NS}(:,neuron);
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
%             kernel = fit.beta(:,end);
%             if s==1 kernels(NN,:)=kernel; end
        end           
    end
    NNstart=NN;
    fprintf('Mouse #%d \n',NS);
end
toc

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

includedModels = 1:length(valueShuffles);
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
b.CData(e,:)=[0.7 0 1]; %trial type
end
b.CData(1,:)=[0 0.7 1]; %value
b.CData(mdlPosition(end),:)=[0.6 0.6 0.6]; %non-specific
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

category=9*ones(length(improvement),1);
category(cueK&~behavior)=4; %cue neurons, other odor coding schemes
category(cueK&~behavior&bestModel==1)=1; %value
category(cueK&~behavior&ismember(bestModel,valRank(2:17)))=2; %value like
category(cueK&~behavior&bestModel==max(bestModel))=3; %untuned

% visual depiction of shuffles
colormap('summer');
subplot(5,1,4);
imagesc(1:length(includedModels),1:6,valueShuffles(includedModels,:)',[0 1]);
xticks([]);
yticks(1:6);
yticklabels({'CS+','CS+','CS50','CS50','CS-','CS-'});
title('Value assigned to each odor in each shuffle');
colorbar;

% visual depiction of shuffles, reordered by correlation
subplot(5,1,5);
includedShuffles = valueShuffles(includedModels,:)';

imagesc(1:length(includedModels),1:6,includedShuffles(:,valRank),[0 1]);
xticks([]);
yticks(1:6);
yticklabels({'CS+','CS+','CS50','CS50','CS-','CS-'});
title('Value assigned to each odor in each shuffle, ordered by similarity to value');
colorbar;

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
thisActivity=dffPSTH3{1,1};
cueResp=mean(thisActivity(:,sortBins),2);
cats={1,2,3};

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
        
        
        activity = dffPSTH3{cue,os}(sel,:);
        activity = activity(sortOrder,:);
        

        
        imagesc(binTimes(plotBins),[1 size(activity,1)],activity(:,plotBins),[-3 3]);%[min(min(cspActivity)) max(max(cspActivity))]);
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
diractivity = mean(cat(3,dffPSTH6{(1-1)*2+2,1}(:,plotBins),dffPSTH6{(1-1)*2+1,2}(:,plotBins)),3);
direction=sign(max(diractivity,[],2)-abs(min(diractivity,[],2)));

    for d=1:2
        subplot(6,3,3+(ct-1)*6+(d-1)*3);
        
        hold on;
        for cue=1:3
            for os=1:2
                
                activity = dffPSTH6{(cue-1)*2+os,os};
                
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
        if d==2
            xlabel('seconds from odor onset');
            ylabel('z-score');
        else
            xticks([]);
        end
        
        %if reg>1 yticks([]); end
        plot([0 0],[-1 5],'color','k','linewidth',0.75);
        plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
        if d==1 & ct==1 axis([plotWindow -0.25 5]); end
        if d==1 & ct==2 axis([plotWindow -0.25 5]); end
        if d==2 axis([plotWindow -0.75 0.5]); end
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