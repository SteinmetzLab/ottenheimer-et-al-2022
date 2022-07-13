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

plotrange = [-1.5 6];
%smoothing filter for licking PSTH
licksmoothbins=10; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',smoothing);
lickfilterweights=pdf(halfnormal,0:licksmoothbins);

anticipatoryLicks={};
lickPSTH3={};
lickPSTH4={};
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
        
        %smooth lick traces
        smoothedLicks=[];
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
        

        
    end
    
    %plot PSTHs
    subplot(1,length(sessions),session)
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
    
    legend([p{:}],'CS+','CS50(r)','CS50(u)','CS-','location','NW');
    title(['Odor set ' num2str(session)]);
    
    
end

figure;
for session=1:length(sessions)
    hold on;
    errorbar([1 2],mean(anticipatoryLicks{1,1}),nanste(anticipatoryLicks{1,1},1),'linewidth',1.5,'color',[0.1 0.6 0.2]);
    errorbar([1 2],mean(anticipatoryLicks{2,1}),nanste(anticipatoryLicks{2,1},1),'linewidth',1.5,'color',[0.4 0.1 0.4]);
    errorbar([1 2],mean(anticipatoryLicks{3,1}),nanste(anticipatoryLicks{3,1},1),'linewidth',1.5,'color',[0.3 0.3 0.3]);
    plot([1 2] + (session-1)*2.2,anticipatoryLicks{1,1},'color',[0.1 0.6 0.2],'linewidth',0.5); 
    plot([1 2] + (session-1)*2.2,anticipatoryLicks{2,1},'color',[0.4 0.1 0.4],'linewidth',0.5); 
    plot([1 2] + (session-1)*2.2,anticipatoryLicks{3,1},'color',[0.3 0.3 0.3],'linewidth',0.5);
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

analysisWindow=[-1 6.5];
binsPerTrial=diff(analysisWindow)/binSize;
predictornames{1}='CS+';
predictornames{2}='CS50';
predictornames{3}='CS-';
predictornames{4}='licks';
predictornames{5}='boutStart';
predictornames{6}='lickreward';
predictornames{7}='lickRate';

%for shifting lick vector in time
%each shift is binSize*2
numShiftsB=3; %making licks happen later than reality
numShiftsF=2; %making licks happen earlier than reality

windows={};
windows{1}=[0 5];
windows{2}=[0 5];
windows{3}=[0 5];
windows{4}=[-0.3 0.3];
windows{5}=[-0.3 2];
windows{6}=[0 analysisWindow(end)-2.5];
windows{7}=[-numShiftsB*binSize (numShiftsF+1)*binSize]; %there are numShifts*2+1 total for this, so fudging it a bit

binPred=[];
binsofar=0;
for win=1:length(windows)
    winbins=diff(windows{win})/binSize;
    binPred(binsofar+1:binsofar+winbins)=win;
    binsofar=binsofar+winbins;
end
lickPreds=ismember(binPred,[4 5 7]);
cuePreds=ismember(binPred,[1 2 3]);
rewPreds=ismember(binPred,[6]);

submodelsels{1}=cuePreds==0;
submodelsels{2}=lickPreds==0;
submodelsels{3}=rewPreds==0;

Xall=cell(length(mice)*length(sessions),1);
Yall=cell(length(mice)*length(sessions),1);
YallSh=cell(length(mice)*length(sessions),1);

dffPSTH3={};
dffPSTH6={};
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
        continuousPredictors=vector1;
        discreteUnfixed.times={cue1,cue2,cue3,lick,boutStart,lickreward};
        discreteUnfixed.windows=windows(1:6);
        discreteUnfixed.binTimes=binTimesGLMincluded(:);
        A=makeA(discreteUnfixed,continuousPredictors);
        Xall{NS,1}=A;
        
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
            
            
            
            %%% prepare for GLM %%%%
            includedActivity=binneddFF(:,includedBins)';
            shOrder=randperm(size(binneddFF,1));
            includedActivitySh=binneddFF(shOrder,includedBins)';
            activityVector=includedActivity(:);
            activityVectorSh=includedActivitySh(:);
            Yall{NS,1}(:,roi)=activityVector;
            YallSh{NS,1}(:,roi)=activityVectorSh;
            
        end
        fprintf('Session #%d \n',NS);
    end
end
toc

%% roi plotting plot PSTHs
map=redblue(256);
figure;
colormap(map);
%heatmaps
cspActivity=dffPSTH3{1,1}; %get the firing rates of neurons of interest
cueWindow=[0 1.5];
cueBins=binTimes>=cueWindow(1) & binTimes<=cueWindow(2);
cueResp=mean(cspActivity(:,cueBins),2);
[xx,sortOrder]=sort(cueResp);
plotWindow=[-0.5 6];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);

colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];

titles={'CS+       ','CS50      ','CS-     '};
sessions={'o1d3','o2d3'};
pn=0;
for cue=1:3
    for session=1:length(sessions)
        sn=sessions{session};

        activity = dffPSTH3{cue,session};
        activity=activity(sortOrder,:);
        pn=pn+1;
        subplot(1,6,pn);
        hold on;
        
        imagesc(binTimes(plotBins),[1 length(cspActivity)],activity(:,plotBins),[-3 3]);%[min(min(cspActivity)) max(max(cspActivity))]);
        ylim([0.5 length(cspActivity)+0.5])
        plot([0 0],[0.5 length(cspActivity)+0.5],'color',colors{cue},'linewidth',0.75);
        plot([2.5 2.5],[0.5 length(cspActivity)+0.5],':','color','k','linewidth',0.25);
        if pn>1
            set(gca,'ytick',[]);
        end
        
        %colorbar;
        title(titles{cue});
        xlabel('seconds from cue');
    end
    
    
    
end


%% get reduced rank predictors

[~, bR, R2] = CanonCor2all(Yall, Xall);

%% perform regression with reduced rank predictors
components=20;
folds=4;
NS=0;
varExp=NaN(length(roiMouse),length(submodelsels)+1);
predF=cell(length(roiMouse),length(submodelsels)+1);
kernels=NaN(length(roiMouse),length(R2));
lambdaVal=NaN(length(roiMouse),1);
lambda=[0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.5];
NN=0;
for session=1:length(Xall)
    
    A = Xall{session};
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
    
    for roi=1:size(Yall{session},2)
        NN=NN+1;
        y=Yall{session}(:,roi); 
        
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
    fprintf('Session #%d \n',session);
    
end

%% GLM results analysis
improvCutoff=0.02;
overallCutoff=0.02;
sessions={'o1d3','o2d3'};
category={};
sessel{1}=cat(1,true(totalNeurons,1),false(totalNeurons,1));
sessel{2}=cat(1,false(totalNeurons,1),true(totalNeurons,1));

for session=1:2
    varExpSes=varExp(sessel{session},:);
    improvement=-(varExpSes-varExpSes(:,1));

    %get proportions
    category{session}=8*ones(length(improvement),1);
    predictive = varExpSes(:,1)>overallCutoff; %more than 5% of variance explained
    behavior = improvement(:,3)>improvCutoff; %num licks or lick rate
    rewardK = improvement(:,4)>improvCutoff; %static or changing cue
    cueK = improvement(:,2)>improvCutoff; %static or changing cue

    %category{session}(predictive,1)=1;
    category{session}(predictive & behavior,1)=3;
    category{session}(predictive & cueK,1)=1;
    category{session}(predictive & rewardK,1)=5;
    category{session}(predictive & rewardK & cueK,1)=6;
    category{session}(predictive & behavior & cueK,1)=2;
    category{session}(predictive & behavior & rewardK,1)=4;
    category{session}(predictive & behavior & cueK & rewardK,1)=7;
    
 
    %mouse
    mice=unique(roiMouse);
    for mouse=1:length(mice)
        mouseSel=strcmp(roiMouse,mice(mouse));

        varianceMouse(mouse,session)=mean(varExpSes(mouseSel,1));
        for variable=1:3
            uniqueVarMouse{variable}(mouse,session)=mean(varExpSes(mouseSel,1)-varExpSes(mouseSel,1+variable));
        end
    end
    

end


figure('position',[100 100 1200 900]);
subplot(3,4,1);
hold on;

expLim=[0 0.2];
aloLim=[-0.01 0.06];
errorbar(1:length(sessions),mean(varianceMouse),nanste(varianceMouse,1),'color',[0 0 0],'linewidth',1);
plot(1:length(sessions),varianceMouse,'color',[0.2 0.2 0.2]);
p=anova1(varianceMouse,[],'off');
if p<0.05 text(2,0.9*expLim(2),'*','fontsize',20); end
xticks(1:3);
xticklabels({'1','2','3'});
xlim([0.5 3.5]);
ylim(expLim);
ylabel('mean R^2');
xlabel('day');

subplot(3,4,2);
vnames={'cues','licks','reward'};
vcolors={[0 0.4 0.9],[0.9 0.2 0],[0.1 0.7 0]};
for variable=1:3
    subplot(3,4,2);
    hold on;
    errorbar((variable-1)*2.5+[1 2],nanmean(uniqueVarMouse{variable}),nanste(uniqueVarMouse{variable},1),'color',vcolors{variable},'linewidth',1);
    plot((variable-1)*2.5+[1 2],uniqueVarMouse{variable},'color',vcolors{variable},'linewidth',1);

end
xlim([0.5 7.5]);
plot([0.5 7.5],[0 0],':','color','k','linewidth',0.75);
ylim(aloLim);
ylabel('unique R^2');
xticks(1:3)
xticklabels(vnames)
xtickangle(45)

subplot(3,4,3);
barData=zeros(8,2);
for session=1:length(sessions)
    for ct=1:8
        barData(ct,session)=sum(category{session}==ct)/length(category{session});
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
xlim([0 6]);
legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
xlabel('day');

%% glm performance across tracked cells

%take kernels from session 1 and test on both sessions
pn=0;
sessions={'o1d3','o2d3'};
%make the matrix for each neuron
varExpTr=NaN(totalNeurons,2);
varExpTrSh=NaN(totalNeurons,2);
neuronnum=[1:totalNeurons]';
for neuron=1:totalNeurons
    mouseNum=find(strcmp(mice,roiMouse(neuron)));
    mouse=mice{mouseNum};
    roi=find(neuronnum(strcmp(roiMouse,mouse))==neuron);
    ks=kernels(neuron,:)';
    testC=0;
    for sessionTest=1:length(sessions)
        NST=(sessionTest-1)*length(mice)+mouseNum;
        AT=Xall{NST};
        y=Yall{NST}(:,roi);
        ySh=YallSh{NST}(:,roi);
        
        predSp = AT*ks;
        varExpTr(neuron,sessionTest) = corr(y,predSp);
        varExpTrSh(neuron,sessionTest) = corr(ySh,predSp);
        
        
    end
end

%% plot GLM across days

%average correlation with other days



meanCorr=[];
meanCorrSh=[];
for mouse=1:length(mice)
    mouseSel=strcmp(roiMouse,mice{mouse});
    meanCorr(mouse,:)=nanmean(varExpTr(mouseSel,:));
    meanCorrSh(mouse,:)=nanmean(varExpTrSh(mouseSel,:));
end


subplot(2,2,2);
hold on;

errorbar(1:2,mean(meanCorr),nanste(meanCorr,1),'linewidth',1.5,'color','k');
errorbar(1:2,mean(meanCorrSh),nanste(meanCorr,1),'linewidth',1.5,'color',[0.5 0.5 0.5]);
plot(1:2,meanCorr,'linewidth',1,'color',[0.2 0.2 0.2]);
plot(1:2,meanCorrSh,'linewidth',1,'color',[0.8 0.8 0.8]);

axis([0 3 0 0.5]);
ylabel('true~predicted corr');
xlabel('day');
legend('true','shuffle');
%title('train late, test early');

data=cat(1,meanCorr,meanCorrSh);
ml=repmat([1;2;3;4;5],2,2);
dl=repmat([1:2],10,1);
gl=cat(1,ones(5,2),2*ones(5,2));
% [p,tbl,stats]=anovan(data(:),{dl(:),gl(:)},'varnames',{'section','shuffle'},'model','interaction');
% c=multcompare(stats,'dimension',[1 2],'ctype','bonferroni');
%if pv(1)<0.05 text(5,0.38,'*','fontsize',20); end

% unique variance across tracked cells

%sort based off of category on day 3
catBase=category{1};

aloLim=[-0.02 0.12];
cats=[1 3 2 4];
catcolors={[0 0.4 0.9];... %cue
[0.9 0.2 0];... %lick
[0.6 0 0.6];... %both
[0.7 0.7 0.7]};
for ct=1:length(cats)
    uniqueVarMouseTracked={};
entry=0;
for session=1:length(sessions)
    varExpSes=varExp(sessel{session},:);
    sn=sessions{session};
    mice=unique(roiMouse);
    
    entry=entry+1;
    for variable=1:3
        
        for mouse=1:length(mice)
            mouseSel=strcmp(roiMouse,mice(mouse))&catBase==cats(ct);
            if cats(ct)==4 mouseSel=strcmp(roiMouse,mice(mouse))&catBase>3; end
            uniqueVarMouseTracked{variable}(mouse,entry)=mean(varExpSes(mouseSel,1)-varExpSes(mouseSel,1+variable));
        end
    end
end





vnames={'cues','licks','reward'};
vcolors={[0 0.4 0.9],[0.9 0.2 0],[0.1 0.7 0]};
for variable=1:3
    
    subplot(2,1,2);
    hold on;
        entries=1:2;
        xvals=entries+(ct-1)*2.5;
        errorbar(xvals,mean(uniqueVarMouseTracked{variable}(:,entries)),nanste(uniqueVarMouseTracked{variable}(:,entries),1),'color',vcolors{variable},'linewidth',1);

end

xticks(1:2);
xticklabels({'set 1','set 2'});
xtickangle(45);
xlim([0 xvals(end)+1]);
ylim(aloLim);
if variable==1 ylabel('unique R^2'); end
title(vnames(variable));
plot([0.5 xvals(end)+0.5],[0 0],':','color','k');

data=cat(1,uniqueVarMouseTracked{:});
ml=repmat([1;2;3;4;5],3,2);
sl=repmat([1:2],15,1);
gl=cat(1,ones(5,2),2*ones(5,2),3*ones(5,2));
% [p,tbl,stats]=anovan(data(:),{sl(:),gl(:)},'varnames',{'section','variable'},'model','interaction');
% c=multcompare(stats,'dimension',[1 2],'ctype','bonferroni');



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
%% nanste function
%standard error, omitting NaNs
function  ste=nanste(dat,dimension)
ste=nanstd(dat,[],dimension)./sqrt(sum(~isnan(dat),dimension));
end