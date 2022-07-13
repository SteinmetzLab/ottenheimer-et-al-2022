% Script to analyze olfactory conditioning data across days 1-3
githubDir = 'D:\GitHub';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'steinmetz-et-al-2019')))
direc = 'D:\GitHub\ottenheimer-et-al-2022\Imaging';
addpath(direc);

%% behavior summary
mice={'pl01','pl02','pl03','pl08','pl10','pl11','pl15','pl16'};
sessions={'o1d1','o1d2','o1d3'};

%parameters
binSize = 0.1;
window = [-3 10];
binEdges = window(1):binSize:window(2);
binTimes = window(1)+binSize/2:binSize:window(2)-binSize/2;

figure;

%parameters
ymax=6.5;
stimDuration = 1.5;
rewOnset = 2.5;
smoothing=8; %controls standard deviation of smoothing filter
window = [-2 10];
plotrange = [-0.5 6];

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

        %get licks/bin
        binnedLicks=NaN(length(cue),length(binEdges)-1);
        antLicks=NaN(length(cue),1);
        preLicks=NaN(length(cue),1);
        for trial = 1:length(cue)
            lickTimesRel = (lick - cue(trial));
            binnedLicks(trial,:)=histcounts(lickTimesRel,binEdges)/binSize;
            antLicks(trial,1) = sum(lickTimesRel>0 & lickTimesRel<2.5);
            preLicks(trial,1) = sum(lickTimesRel>-2.5 & lickTimesRel<0);
        end
       
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
            sel = ismember(round(cue),round(selections{condition})) & trialNo<=60; % TRIAL SPLIT HERE
            lickPSTH3{condition,1}(mouse,:)=nanmean(smoothedLicks(sel,:));
            anticipatoryLicks{condition,session}(mouse,1)=nanmean(antLicks(sel))-nanmean(preLicks(sel));
            
            sel = ismember(round(cue),round(selections{condition})) & trialNo>trialNo(end)-60; % TRIAL SPLIT HERE
            lickPSTH3{condition,2}(mouse,:)=nanmean(smoothedLicks(sel,:)); 
            anticipatoryLicks{condition,session}(mouse,2)=nanmean(antLicks(sel))-nanmean(preLicks(sel));
        end         
        
        %4 conditions
        selections{1,1}=cue1; %cs+
        selections{2,1}=cue2r; %reward+
        selections{3,1}=cue2u; %reward-
        selections{4,1}=cue3; %cs- 
        for condition=1:length(selections)
            sel = ismember(round(cue),round(selections{condition})) & trialNo<=60; % TRIAL SPLIT HERE
            lickPSTH4{condition,1}(mouse,:)=nanmean(smoothedLicks(sel,:));
            sel = ismember(round(cue),round(selections{condition})) & trialNo>trialNo(end)-60; % TRIAL SPLIT HERE
            lickPSTH4{condition,2}(mouse,:)=nanmean(smoothedLicks(sel,:));
        end
        
    end
    
    %plot PSTHs
    for sessionStage=1:2
        subplot(length(sessions),2,sessionStage+2*(session-1))
        hold on

        colors{1,1}=[0.1 0.6 0.2];
        colors{2,1}=[0.4 0.1 0.4];
        colors{3,1}=[0.3 0.3 0.3];
        xvals = find(binTimes>=plotrange(1) & binTimes<=2.5);
        for condition=1:3

            %get values
            psth=nanmean(lickPSTH3{condition,sessionStage});
            sem=nanste(lickPSTH3{condition,sessionStage},1); %calculate standard error of the mean
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
            psth=nanmean(lickPSTH4{condition,sessionStage});
            sem=nanste(lickPSTH4{condition,sessionStage},1); %calculate standard error of the mean
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
        if sessionStage==1 title('early'); end
        if sessionStage==2 title('late'); end
    end
    
end

figure;
for session=1:length(sessions)
    hold on;
    errorbar([1 2] + (session-1)*2.2,mean(anticipatoryLicks{1,session}),nanste(anticipatoryLicks{1,session},1),'linewidth',1.5,'color',[0.1 0.6 0.2]);
    errorbar([1 2] + (session-1)*2.2,mean(anticipatoryLicks{2,session}),nanste(anticipatoryLicks{2,session},1),'linewidth',1.5,'color',[0.4 0.1 0.4]);
    errorbar([1 2] + (session-1)*2.2,mean(anticipatoryLicks{3,session}),nanste(anticipatoryLicks{3,session},1),'linewidth',1.5,'color',[0.3 0.3 0.3]);
    plot([1 2] + (session-1)*2.2,anticipatoryLicks{1,session},'color',[0.1 0.6 0.2],'linewidth',0.5); 
    plot([1 2] + (session-1)*2.2,anticipatoryLicks{2,session},'color',[0.4 0.1 0.4],'linewidth',0.5); 
    plot([1 2] + (session-1)*2.2,anticipatoryLicks{3,session},'color',[0.3 0.3 0.3],'linewidth',0.5);
    ylabel('\Delta anticipatory licks');
    xticks([1 2 3.2 4.2 5.4 6.4]);
    xticklabels({'early','late','early','late','early','late'});
    text(1.2+ (session-1)*2.2,13,['day ' num2str(session)]);
    plot([0 8],[0 0],':','color','k','linewidth',0.75);
    xtickangle(45);
    xlim([0.5 7]);
end

%% get predictors and activity from all sessions and ROIs
mice={'pl01','pl02','pl03','pl08','pl10','pl11','pl15','pl16'};
sessions={'o1d1','o1d2','o1d3'};
neuCoeff=0.7; %coefficient for subtracting neuropil, this is default

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
conPreds=ismember(binPred,8);

submodelsels{1}=cuePreds==0;
submodelsels{2}=lickPreds==0;
submodelsels{3}=rewPreds==0;

Xall=cell(length(mice)*length(sessions),1);
Yall=cell(length(mice)*length(sessions),1);
YallSh=cell(length(mice)*length(sessions),1);
dffPSTH3=struct;
dffPSTH4=struct;
NS=0;
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
        load(fullfile(direc,mn,sn,'Fall.mat'));
        
        % remove doubles and ROIs on edge of FOV
        iscell = processROIs(iscell,stat,F);
        
        numROIs=sum(iscell(:,1));
        indx=find(iscell(:,1));

        for roi=1:numROIs
            
            NN=NN+1; 
            roiNum=indx(roi);
%             correctedF=F(roiNum,:)-(neuCoeff * Fneu(roiNum,:));
%             medianF=medfilt1(double(correctedF),7);        

            
            
            decSpks=spks(roiNum,:);
            
            %smooth with the same filter used for ephys
            spikeRateSmooth=NaN(1,length(decSpks));
            for l=1:length(spikeRateSmooth)
                spikeRateSmooth(1,l)=sum(decSpks(1,l-min([l-1 smoothbinsTrl]):l).*fliplr(filterweightsTrl(1:min([l smoothbinsTrl+1]))))/sum(filterweightsTrl(1:min([l smoothbinsTrl+1])));
            end            
            
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
            
            %to avoid silly z-scores
            neuronStd(NN,1)=std(binnedActivity(:,binTimes<0),0,'all');
            if std(binnedActivity(:,binTimes<0),0,'all')>=1
                binneddFF=(binnedActivity-mean(binnedActivity(:,binTimes<0),'all')) / std(binnedActivity(:,binTimes<0),0,'all');
            else
                binneddFF=(binnedActivity-mean(binnedActivity(:,binTimes<0),'all'));
            end
            
%             
%             %convert to dF/F
%             binneddFF = (binnedActivity - median(binnedActivity(:,binTimes<0),2)) / median(F(roiNum,:));            

            %3 conditions
            selections={};
            selections{1,1}=cue1; %cs+
            selections{2,1}=cue2; %cs50
            selections{3,1}=cue3; %cs-
            for condition=1:length(selections)
                sel = ismember(round(cue),round(selections{condition})) & trialNo<=60; %% TRIAL CUTOFF HERE
                dffPSTH3.(sn){condition,1}(NN,:)=nanmean(binneddFF(sel,:));

                sel = ismember(round(cue),round(selections{condition})) & trialNo>trialNo(end)-60; %% TRIAL CUTOFF HERE
                dffPSTH3.(sn){condition,2}(NN,:)=nanmean(binneddFF(sel,:)); 
                
                sel = ismember(round(cue),round(selections{condition}));
                dffPSTH3.(sn){condition,3}(NN,:)=nanmean(binneddFF(sel,:)); 
            end                    
            
            %4 conditions
            selections{1,1}=cue1; %cs+
            selections{2,1}=cue2r; %reward+
            selections{3,1}=cue2u; %reward-
            selections{4,1}=cue3; %cs- 
            for condition=1:length(selections)
                sel = ismember(round(cue),round(selections{condition}));
                dffPSTH4.(sn){condition,1}(NN,:)=nanmean(binneddFF(sel,:));                
            end
            
            %6 conditions
            selections={};
            selections{1,1}=cue1; %cs+
            selections{2,1}=cue1; %cs+
            selections{3,1}=cue2; %cs50
            selections{4,1}=cue2; %cs50
            selections{5,1}=cue3; %cs-
            selections{6,1}=cue3; %cs-
            for condition=1:length(selections)
                sel = ismember(round(cue),round(selections{condition})) & rem(trialNo,2)==rem(condition,2);
                dffPSTH6.(sn){condition,1}(NN,:)=nanmean(binneddFF(sel,:));   
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
figure('position',[100 100 500 1200]);
colormap(map);
%heatmaps
cueWindow=[0 1.5];
cueBins=binTimes>=cueWindow(1) & binTimes<=cueWindow(2);


colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];

plotWindow=[-0.5 6];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);
titles={'CS+','CS50','CS-'};
sessions={'o1d1','o1d2','o1d3'};
en=0;
for session=1:length(sessions)
    sn=sessions{session};
    cspActivity=dffPSTH3.(sn){1,2}; %get the firing rates of neurons of interest
    cueResp=mean(cspActivity(:,cueBins),2);
    
    [xx,sortOrder]=sort(cueResp);

    pn=0;
    for sessionStage=1:2
        en=en+1;
        for cue=1:3
            
            activity = dffPSTH3.(sn){cue,sessionStage};
            activityS=activity(sortOrder,:);
            pn=pn+1;
            subplot(length(sessions),6,(session-1)*6+pn);
            hold on;
            
            imagesc(binTimes(plotBins),[1 length(cspActivity)],activityS(:,plotBins),[-3 3]);%[min(min(cspActivity)) max(max(cspActivity))]);
            ylim([0.5 length(cspActivity)+0.5])
            plot([0 0],[0.5 length(cspActivity)+0.5],'color',colors{cue},'linewidth',0.75);
            plot([2.5 2.5],[0.5 length(cspActivity)+0.5],':','color','k','linewidth',0.25);
            
            if pn>1
                set(gca,'ytick',[]);
            else
                yticks([0 length(cspActivity)]);
                ylabel(sprintf('Individual ROIs, day %g',session));
            end
            if session==1 title(titles{cue}); end
            xticks([]);
            
            
        end
        
        %colorbar;
        if session==3
        xlabel('seconds from cue');
        xticks([0 2.5 6]);
        end
        
    end
end

%% get reduced rank predictors

[~, bR, R2] = CanonCor2all(Yall, Xall);

%% perform regression with reduced rank predictors
components=20;
folds=4;
NS=0;
varExp=struct;
roiMouse=struct;
predF=struct;
kernels=struct;
lambdaVal=struct;
lambda=[0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.5];

tic

for session=1:length(sessions)
    sn=sessions{session};
    NN=0;
    for mouse=1:length(mice)
        NS=NS+1;
        mn=mice{mouse};
        
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
        
        for roi=1:size(Yall{NS},2)
            NN=NN+1;
            roiMouse.(sn){NN,1}=mn;
            y=Yall{NS}(:,roi)-mean(Yall{NS}(:,roi));
            
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
            [varExp.(sn)(NN,1),ind] = max(varExpLam);
            thisLambda=lambda(ind);
            lambdaVal.(sn)(NN,1)=thisLambda;

            pred=predLam(:,ind);
            predF.(sn){NN,1}=pred;
  
            %full data to get kernels
            fitK=lassoglm(rA,y,'normal','alpha',0.5,'lambda',thisLambda);
            kernel=bR(:,1:components)*fitK;
            kernels.(sn)(NN,1:length(kernel))=kernel;
            
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
                predF.(sn){NN,1+sub}=pred;
                varExp.(sn)(NN,1+sub) = 1- var(y-pred)/var(y);
                
            end
            
        end
        fprintf('Session #%d \n',NS);
    end
end

toc

%% GLM results analysis
improvCutoff=0.02;
overallCutoff=0.02;
sessions={'o1d1','o1d2','o1d3'};
category={};
varianceMouse=[];
uniqueVarMouse={};
for session=1:length(sessions)
    sn=sessions{session};
    improvement=-(varExp.(sn)-varExp.(sn)(:,1));

    %get proportions
    category{session}=8*ones(length(improvement),1);
    predictive = varExp.(sn)(:,1)>overallCutoff; %more than 5% of variance explained
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
    mice=unique(roiMouse.(sn));
    for mouse=1:length(mice)
        mouseSel=strcmp(roiMouse.(sn),mice(mouse));
        varianceMouse(mouse,session)=mean(varExp.(sn)(mouseSel,1));
        for variable=1:3
            uniqueVarMouse{variable}(mouse,session)=mean(varExp.(sn)(mouseSel,1)-varExp.(sn)(mouseSel,1+variable));
        end
    end

end


figure('position',[100 100 1200 900]);
subplot(3,4,1);
hold on;

expLim=[0 0.2];
aloLim=[-0.01 0.07];
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
barData=[];
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

vnames={'cues','licks','reward'};
vcolors={[0 0.4 0.9],[0.9 0.2 0],[0.1 0.7 0]};

for variable=1:3

    subplot(6,8,13+variable);
    hold on;
    errorbar(1:length(sessions),mean(uniqueVarMouse{variable}),nanste(uniqueVarMouse{variable},1),'color',vcolors{variable},'linewidth',1);
    plot(1:length(sessions),uniqueVarMouse{variable},'color',vcolors{variable}*1.1);
    p=anova1(uniqueVarMouse{variable},[],'off');
    if p<0.05 text(2,0.9*aloLim(2),'*','fontsize',20); end
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0.5 3.5]);
    ylim(aloLim);
    if variable==1 ylabel('unique R^2'); end
    xlabel('day');
  
end

    
%% perform regression with reduced rank predictors on session sections
sections=3;
components=20;
folds=10;
NS=0;
varExpSec=struct;
kernelsSec=struct;
lambda=[0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.5];
tic
for session=1:length(sessions)
    sn=sessions{session};
    NN=0;
    for mouse=1:length(mice)
        NS=NS+1;
        mn=mice{mouse};
        
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
        
        
        for roi=1:size(Yall{NS},2)
            NN=NN+1;
            y=Yall{NS}(:,roi)-mean(Yall{NS}(:,roi));

            for section=1:sections
                secSel=false(length(y),1);
                secSel(1+floor((section-1)*length(A)/sections):floor(section*length(A)/sections))=true;
                

                %cross-validated variance explained
                predLam=NaN(length(secSel),length(lambda));
                for fold=1:folds
                    train=trains{fold}&secSel;
                    test=train==0&secSel;
                    fitK=lassoglm(rA(train,:),y(train),'normal','alpha',0.5,'lambda',lambda);
                    kernel=bR(:,1:components)*fitK;
                    predLam(test,:)=A(test,:)*kernel;
                end
                
                %get best lambda
                for lam=1:length(lambda)
                    varExpLam(1,lam)=1- var(y(secSel)-predLam(secSel,lam))/var(y(secSel));
                end
                [varExpSec.(sn){section}(NN,1),ind] = max(varExpLam);
                thisLambda=lambda(ind);
                
                
                %full data to get kernels
                fitK=lassoglm(rA(secSel,:),y(secSel),'normal','alpha',0.5,'lambda',thisLambda);
                kernel=bR(:,1:components)*fitK;
                kernelsSec.(sn){section}(NN,1:length(kernel))=kernel;
                
                %submodels to find unique variance and total variance for each
                %variable
                for sub=1:length(submodelsels)
                    pred=NaN(length(secSel),1);
                    for fold=1:folds
                        train=trains{fold}&secSel;
                        test=train==0&secSel;
                        fitK=lassoglm(rAs{sub}(train,:),y(train),'normal','alpha',0.5,'lambda',thisLambda);
                        kernel=bR(:,1:components)*fitK;
                        pred(test)=As{sub}(test,:)*kernel;
                    end
                    predAll=pred;
                    pred(isnan(pred))=[];
                    varExpSec.(sn){section}(NN,1+sub) = 1- var(y(secSel)-pred)/var(y(secSel));
                    
                end
            end
        end
        fprintf('Session #%d \n',NS);
    end
end
toc
%% GLM results for sections
plotTraces=true;
improvCutoff=0.02;
overallCutoff=0.02;
sessions={'o1d1','o1d2','o1d3'};
categorys={};

figure('position',[100 100 1200 900]);
entry=0;
cuesMouse=[];
rewMouse=[];
lickMouse=[];
varianceMouse=[];
uniqueVarMouse={};

for session=1:length(sessions)
    sn=sessions{session};
    for section=1:3
        entry=entry+1;
        improvement=-(varExpSec.(sn){section}-varExpSec.(sn){section}(:,1));

        %get proportions
        predictive = varExpSec.(sn){section}(:,1)>overallCutoff; %more than 5% of variance explained
        behavior = improvement(:,3)>improvCutoff; %num licks or lick rate
        rewardK = improvement(:,4)>improvCutoff; %static or changing cue
        cueK = improvement(:,2)>improvCutoff; %static or changing cue
        
        %category{session}(predictive,1)=1;
        categorys{session,section}=8*ones(length(improvement),1);        
        categorys{session,section}(predictive & behavior,1)=3;
        categorys{session,section}(predictive & cueK,1)=1;
        categorys{session,section}(predictive & rewardK,1)=5;
        categorys{session,section}(predictive & rewardK & cueK,1)=6;
        categorys{session,section}(predictive & behavior & cueK,1)=2;
        categorys{session,section}(predictive & behavior & rewardK,1)=4;
        categorys{session,section}(predictive & behavior & cueK & rewardK,1)=7;
        
        
        %mouse
        mice=unique(roiMouse.(sn));
        for mouse=1:length(mice)
            mouseSel=strcmp(roiMouse.(sn),mice(mouse));
            cuesMouse(mouse,entry)=sum(predictive&cueK&mouseSel)/sum(mouseSel);
            rewMouse(mouse,entry)=sum(predictive&rewardK&mouseSel)/sum(mouseSel);
            lickMouse(mouse,entry)=sum(predictive&behavior&mouseSel)/sum(mouseSel);
            varianceMouse(mouse,entry)=mean(varExpSec.(sn){section}(mouseSel,1));
            for variable=1:3
                uniqueVarMouse{variable}(mouse,entry)=mean(varExpSec.(sn){section}(mouseSel,1)-varExpSec.(sn){section}(mouseSel,1+variable));
            end
        end
   
    end
    
end
    
subplot(3,4,1);
hold on;

expLim=[-0.025 0.15];
aloLim=[-0.01 0.06];
for d=1:3
    entries=(d-1)*3+1:d*3;
    errorbar(entries,mean(varianceMouse(:,entries)),nanste(varianceMouse(:,entries),1),'color',[0 0 0],'linewidth',1);
    plot(entries,varianceMouse(:,entries),'color',[0.2 0.2 0.2]);
    [p,tbl]=anova1(varianceMouse(:,entries),[],'off');
    if p<0.05 text(entries(2),0.5*expLim(2),'*','fontsize',20); end
end
xticks(1:3);
xticklabels({'early','middle','late'});
xtickangle(45);
xlim([0.5 9.5]);
plot([0.5 9.5],[0 0],':','color','k');
ylim(expLim);
ylabel('mean R^2');


subplot(3,4,2);
barData=[];
session=1;
for section=1:3
    for ct=1:8
        barData(ct,section)=sum(categorys{session,section}==ct)/length(categorys{session,section});
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
xlabel('section');

vnames={'cues','licks','reward'};
vcolors={[0 0.4 0.9],[0.9 0.2 0],[0.1 0.7 0]};
subplot(3,2,2);
hold on;
for variable=1:3
    

    for d=1:3
        entries=(d-1)*3+1:d*3;
        xvals=entries+(variable-1)*10;
    errorbar(xvals,mean(uniqueVarMouse{variable}(:,entries)),nanste(uniqueVarMouse{variable}(:,entries),1),'color',vcolors{variable},'linewidth',1);
    plot(xvals,uniqueVarMouse{variable}(:,entries),'color',vcolors{variable}*1.1);
    [p,tbl,stat]=anova1(uniqueVarMouse{variable}(:,entries),[],'off');
    %multcompare(stat);
    if p<0.05 text(xvals(2),0.5*aloLim(2),'*','fontsize',20); end
    
    end
    


end

    xticks(1:3);
    xticklabels({'early','middle','late'});
    xtickangle(45);
    xlim([0.5 xvals(end)+0.5]);
    ylim(aloLim);
    if variable==1 ylabel('unique R^2'); end
    title(vnames(variable));
    plot([0.5 xvals(end)+0.5],[0 0],':','color','k');
%% plot traces of cue and lick cells


colors{1,1}=[0.1 0.6 0.2];
colors{2,1}=[0.4 0.1 0.4];
colors{3,1}=[0.3 0.3 0.3];

plotWindow=[-1 6];
plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);
xvals=binTimes(plotBins);

cats=[1 3 2 4];
catnames={'cue cells','lick cells','both cells','others'};
session=1;
sn=sessions{session};

catcolors={[0 0.4 0.9];... %cue
[0.9 0.2 0];... %lick
[0.6 0 0.6];... %both
[0.7 0.7 0.7]};

figure;

directions=[1 -1];
for condition=1:3
    diractivity = dffPSTH3.(sn){condition,2};%(:,cueBins);
    direction=sign(max(diractivity,[],2)-abs(min(diractivity,[],2)));
    
    for d=1:2
        
        
        
        for el=1:2
            subplot(4,6,condition+(d-1)*6+(el-1)*3);%+(reg-1)*6);
            hold on;
            
            activity = dffPSTH3.(sn){condition,el};
            
            for reg=1:length(cats)
                regsel=ismember(categorys{session,3},cats(reg));
                if cats(reg)==4 regsel=categorys{session,3}>3; end
                
                %get values
                psth=nanmean(activity(regsel&direction==directions(d),plotBins));
                sem=nanste(activity(regsel&direction==directions(d),plotBins),1); %calculate standard error of the mean
                up=psth+sem;
                down=psth-sem;
                
                %plotting
                plot(binTimes(plotBins),psth,'Color',catcolors{reg,1},'linewidth',0.75);
                patch([xvals,xvals(end:-1:1)],[up,down(end:-1:1)],catcolors{reg,1},'EdgeColor','none');alpha(0.2);
                
            end
            if d==1 axis([plotWindow -0.25 1.75]); end
            if d==2 axis([plotWindow -0.75 0.5]); end
            plot([2.5 2.5],[-2 4],'color','k','linewidth',0.75);
            patch([0 0 1.5 1.5],[-2 4 4 -2],colors{condition},'edgecolor','none');alpha(0.3);
            plot(plotWindow,[0 0],':','color','k','linewidth',0.5);
            
            if condition==1 & d==2
                xlabel('seconds from odor onset');
                ylabel('z-score');
            else
                xticks([]);
            end
            
            if condition>1 yticks([]); end
        end
        
    end
    
    
    
end



%% load in manually tracked cells for days 1-3


trackingDir='D:\Dropbox\Manuscripts\PL\ottenheimerEtAl2022\Imaging\ROIs\Masks';
mice={'pl01','pl02','pl03','pl08','pl10','pl11','pl15','pl16'};
sessions={'o1d1','o1d2','o1d3'};
trackedROIids={};
for session=1:length(sessions)
    sn=sessions{session};
    NN=0;
    trackID=0;
    roiID=[];
    for mouse=1:length(mice)
        mn=mice{mouse};
        
        %get ROI data
        load(fullfile(direc,mn,sn,'Fall.mat'));      
        iscell=processROIs(iscell,stat,F);
        
        roiID=cat(1,roiID,zeros(sum(iscell(:,1)),1));
        %get location of each tracked ROI for this session
        sessionName=append(mn,sn);
        trackCoords=readNPY(fullfile(trackingDir,append(sessionName,'ROImasks.npy')));
        %fprintf('Session %d, Mouse %d, %d ROIs \n',session,mouse,length(trackCoords));
        idx=find(iscell(:,1));
        roix=NaN(sum(iscell(:,1)),1);
        roiy=NaN(sum(iscell(:,1)),1);
        for roi=1:sum(iscell(:,1))
            roix(roi,1)=stat{1,idx(roi)}.med(2);
            roiy(roi,1)=stat{1,idx(roi)}.med(1);
        end
        
        distances=NaN(sum(iscell(:,1)),1);
        trackIDs=[];
        distance=[];
        closestROI=[];
        distance2={};
        closestROI2={};
        for tr=1:length(trackCoords)
            trackID=trackID+1;
            trackIDs(tr,1)=trackID;
            distances=sqrt((roix-trackCoords(1,tr)).^2+(roiy-trackCoords(2,tr)).^2);
            [distance(tr,1),closestROI(tr,1)]=min(distances);
            %get second closest
            [distance2{tr,1},closestROI2{tr,1}]=sort(distances);
            
        end
        instances=histcounts(closestROI,0.5:1:sum(iscell(:,1))+0.5);
        fprintf('%d assigned more than once \n',sum(instances>1));
        while sum(instances>1)>0
            doubles=find(instances>1);
            for dbl=1:length(doubles)
                badInd=find(closestROI==doubles(dbl));
                doubDistances=distance(badInd);
                farAways=find(doubDistances>min(doubDistances));
                if sum(farAways)<1
                    doubDistances2=[];
                    for d=1:length(badInd)
                    doubDistances2(d,1)=distance2{badInd(d)}(2);
                    end
                    [x,farAways]=min(doubDistances2);
                end
                for change=1:length(farAways)
                    closestROI(badInd(farAways(change)))=closestROI2{badInd(farAways(change))}(1);
                    distance(badInd(farAways(change)))=distance2{badInd(farAways(change))}(1);
                    closestROI2{badInd(farAways(change))}(1)=[];
                    distance2{badInd(farAways(change))}(1)=[];
                end
            end
            instances=histcounts(closestROI,0.5:1:sum(iscell(:,1))+0.5);
        end
                
        roiID(NN+closestROI,1)=trackIDs;
        NN=NN+sum(iscell(:,1));
    end
    fprintf('%d tracked from session %d \n',trackID,session);
    trackedROIids{session}=roiID;
    
    
end

%% plot cells, same across all days
plotCD=false;

map=redblue(256);
figure;
colormap(map);

%sorting
cspActivity=dffPSTH3.o1d3{1,2}(trackedROIids{3}>0,:); %get the firing rates of neurons of interest
[~,idOrder]=sort(trackedROIids{3}(trackedROIids{3}>0));
cspActivity=cspActivity(idOrder,:);
cueWindow=[0 1.5];
cueBins=binTimes>=cueWindow(1) & binTimes<=cueWindow(2);
cueResp=mean(cspActivity(:,cueBins),2);
[~,sortOrder]=sort(cueResp);

earlyLate={'early','late'};
sessions={'o1d1','o1d2','o1d3'};
pn=0;
for cue=1:3
    for session=1:length(sessions)
        sn=sessions{session};
        
        plotWindow=[-0.5 6];
        plotBins=binTimes>=plotWindow(1) & binTimes<=plotWindow(2);
        
        titles={'CS+','CS50','CS-'};
        
        for sessionStage=3
            activity = dffPSTH3.(sn){cue,sessionStage}(trackedROIids{session}>0,:);
            [~,idOrder]=sort(trackedROIids{session}(trackedROIids{session}>0));
            activity=activity(idOrder,:);
            activity=activity(sortOrder,:);
            pn=pn+1;
            
            subplot(1,9,pn);
            hold on;
            
            imagesc(binTimes(plotBins),[1 length(cspActivity)],activity(:,plotBins),[-3 3]);%[min(min(cspActivity)) max(max(cspActivity))]);
            ylim([1 length(cspActivity)])
            xlim(plotWindow);
            title(sprintf(append(titles{cue},', day %d      '), session));
            plot([0 0],[0.5 length(sortOrder)+0.5],'color',colors{cue},'linewidth',0.75);
            plot([2.5 2.5],[0.5 length(sortOrder)+0.5],':','color','k','linewidth',0.25);
            
            if pn>1
                set(gca,'ytick',[]);
            end
            
        end
                
        if session==1 xlabel('seconds from odor onset'); end
        
    end
end


%% show one cell's traces on every trial across days

roiIDs=[sortOrder([354 367 359 363])];
plotWindow=[-1 6];
plotBins=find(binTimes>=plotWindow(1) & binTimes<=plotWindow(2));
limits=[-3 3];
for roi=1:length(roiIDs)
    roiID=roiIDs(roi);
    sessions={'o1d1','o1d2','o1d3'};
    mn=roiMouse.o1d1{ismember(trackedROIids{1},roiID)};
    allTrials=[];
    sessBreak=[];
    for session=1:length(sessions)
        sn=sessions{session};
        sessionName=append(mn,sn);
        load(fullfile(direc,mn,sn,append(sessionName,'events.mat')));
        load(fullfile(direc,mn,sn,'Fall.mat'));
        iscell = processROIs(iscell,stat,F);
        numROIs=sum(iscell(:,1));
        indx=find(iscell(:,1));
        sessionROIid=find(trackedROIids{session}(strcmp(roiMouse.(sn),mn),1)==roiID);
        roiNum=indx(sessionROIid);
        %         correctedF=F(roiNum,:)-(neuCoeff * Fneu(roiNum,:));
        %         medianF=medfilt1(double(correctedF),7);
        decSpks=spks(roiNum,:);
        spikeRateSmooth=NaN(1,length(decSpks));
        
        spikeRateSmooth=NaN(1,length(decSpks));
        for l=1:length(spikeRateSmooth)
            spikeRateSmooth(1,l)=sum(decSpks(1,l-min([l-1 smoothbinsTrl]):l).*fliplr(filterweightsTrl(1:min([l smoothbinsTrl+1]))))/sum(filterweightsTrl(1:min([l smoothbinsTrl+1])));
        end
        
        binnedActivity=NaN(length(cue1),length(binEdges)-1);
        numBins=length(binTimes);
        for trial = 1:length(cue1)
            frameTimesRel = (frameTimes - cue1(trial));
            for bin=1:numBins
                binnedActivity(trial,bin)=mean(spikeRateSmooth(frameTimesRel>=binEdges(bin)&frameTimesRel<binEdges(bin+1)));
            end
        end
        %binneddFF = (binnedActivity - mean(binnedActivity(:,binTimes<0),2)) / median(F(roiNum,:)  );
        %to avoid silly z-scores
        if std(binnedActivity(:,binTimes<0),0,'all')>=1
            binneddFF=(binnedActivity-mean(binnedActivity(:,binTimes<0),'all')) / std(binnedActivity(:,binTimes<0),0,'all');
        else
            binneddFF=(binnedActivity-mean(binnedActivity(:,binTimes<0),'all'));
        end
        %binneddFF=binnedActivity;
        allTrials=cat(1,allTrials,binneddFF);
        sessBreak(session,1)=size(allTrials,1);
    end

    %limits=[-max(allTrials,[],'all') max(allTrials,[],'all')];
    subplot(1,length(roiIDs),roi);
    hold on
    colormap(map);
    imagesc(binTimes(plotBins),[1 size(allTrials,1)],allTrials(:,plotBins),limits);%[min(min(cspActivity)) max(max(cspActivity))]);
    for line=1:length(sessBreak)-1
        plot([binTimes(plotBins(1)) binTimes(plotBins(end))],[sessBreak(line)+0.5 sessBreak(line)+0.5],'linewidth',1,'color','k');
    end
    
    set(gca,'ydir','reverse');
    xlim(plotWindow);
    xticks(0:2:6);
    ylim([0.5 size(allTrials,1)+0.5]);
    plot([0 0],[0.5 size(allTrials,1)+0.5],'color',colors{1},'linewidth',0.75);
    plot([2.5 2.5],[0.5 size(allTrials,1)+0.5],':','color','k','linewidth',0.25);
    yticks([]);
end
xlabel('seconds from odor onset');
%% glm performance across tracked cells

%starter model
session=3;


pn=0;
sessions={'o1d1','o1d2','o1d3'};
sections=3;
%make the matrix for each neuron
varExpTr=NaN(max(trackedROIids{1}),9);
varExpTrSh=NaN(max(trackedROIids{1}),9);
trackedMouse={};
for neuron=1:max(trackedROIids{1})
    mouseNum=find(strcmp(mice,roiMouse.o1d1(trackedROIids{1,1}==neuron)));
    NN=find(trackedROIids{session}==neuron);
    mouse=mice{mouseNum};
    trackedMouse{neuron,1}=mouse;
    ks=kernels.(sn)(NN,:)';
    testC=0;
    for sessionTest=1:length(sessions)
        snT=sessions{sessionTest};
        roi=find(trackedROIids{sessionTest}(strcmp(roiMouse.(snT),mouse))==neuron);
        NST=(sessionTest-1)*length(mice)+mouseNum;
        AT=Xall{NST};
        y=Yall{NST}(:,roi);
        ySh=YallSh{NST}(:,roi);
        for sectionTest=1:3
            secSel=false(length(y),1);
            secSel(1+floor((sectionTest-1)*length(AT)/sections):floor(sectionTest*length(AT)/sections))=true;
            testC=testC+1;
            Asec=AT(secSel,:);
            ySec=y(secSel);
            ySecSh=ySh(secSel);
            predSpSec = Asec*ks;
            predSpSecSh = Asec*ks;
            varExpTr(neuron,testC) = corr(ySec,predSpSec);
            varExpTrSh(neuron,testC) = corr(ySecSh,predSpSec);
            
        end
    end
end

%% plot GLM across days

subplot(2,3,1);
barData=[];
for session=3
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
xlim([2 8]);
legend('Cues','Cues,Licks','Licks','Licks,Reward','Reward','Cues,Reward','Cues,Licks,Reward');
xlabel('day');


%average correlation with other days
meanCorr=[];
meanCorrSh=[];
for mouse=1:length(mice)
    mouseSel=strcmp(trackedMouse,mice{mouse});%&catTrack==1;
    meanCorr(mouse,:)=nanmean(varExpTr(mouseSel,:));
    meanCorrSh(mouse,:)=nanmean(varExpTrSh(mouseSel,:));
end


subplot(2,3,3);
hold on;
for day=0:2

errorbar(day*3+1:day*3+3,mean(meanCorr(:,day*3+1:day*3+3)),nanste(meanCorr(:,day*3+1:day*3+3),1),'linewidth',1.5,'color','k');
errorbar(day*3+1:day*3+3,mean(meanCorrSh(:,day*3+1:day*3+3)),nanste(meanCorr(:,day*3+1:day*3+3),1),'linewidth',1.5,'color',[0.5 0.5 0.5]);
plot(day*3+1:day*3+3,meanCorr(:,day*3+1:day*3+3),'linewidth',1,'color',[0.2 0.2 0.2]);
plot(day*3+1:day*3+3,meanCorrSh(:,day*3+1:day*3+3),'linewidth',1,'color',[0.8 0.8 0.8]);
end
axis([0.5 9.5 0 0.5]);
ylabel('true~predicted corr');
xlabel('day');
legend('true','shuffle');
%title('train late, test early');

data=cat(1,meanCorr,meanCorrSh);
ml=repmat([1;2;3;4;5;6;7;8],2,9);
dl=repmat([1:9],16,1);
gl=cat(1,ones(8,9),2*ones(8,9));
% [p,tbl,stats]=anovan(data(:),{dl(:),gl(:)},'varnames',{'section','shuffle'},'model','interaction');
% c=multcompare(stats,'dimension',[1 2],'ctype','bonferroni');
%if pv(1)<0.05 text(5,0.38,'*','fontsize',20); end

% unique variance across tracked cells

%sort based off of category on day 3
catInc=category{3}(trackedROIids{3}>0);
[~,idOrder]=sort(trackedROIids{3}(trackedROIids{3}>0));
catTrack=catInc(idOrder);


aloLim=[-0.02 0.13];
cats=[1 3 2 4];
catcolors={[0 0.4 0.9];... %cue
[0.9 0.2 0];... %lick
[0.6 0 0.6];... %both
[0.7 0.7 0.7]};
for ct=1:length(cats)
    uniqueVarMouseTracked={};
entry=0;
for session=1:length(sessions)
    sn=sessions{session};
    mice=unique(roiMouse.(sn));
    for section=1:3
        entry=entry+1;
        varExpTracked=varExpSec.(sn){section}(trackedROIids{session}>0,:);
        [~,idOrder]=sort(trackedROIids{session}(trackedROIids{session}>0));
        varExpTracked=varExpTracked(idOrder,:);
        for variable=1:3
            
            for mouse=1:length(mice)
                mouseSel=strcmp(trackedMouse,mice(mouse))&catTrack==cats(ct);
                if cats(ct)==4 mouseSel=strcmp(trackedMouse,mice(mouse))&catTrack>3; end
                uniqueVarMouseTracked{variable}(mouse,entry)=mean(varExpTracked(mouseSel,1)-varExpTracked(mouseSel,1+variable));
            end
        end
    end
end

vnames={'cues','licks','reward'};
vcolors={[0 0.4 0.9],[0.9 0.2 0],[0.1 0.7 0]};
for variable=1:3
    
    subplot(2,1,2);
    hold on;
    for d=1:3
        entries=(d-1)*3+1:d*3;
        xvals=entries+(ct-1)*10;
        errorbar(xvals,mean(uniqueVarMouseTracked{variable}(:,entries)),nanste(uniqueVarMouseTracked{variable}(:,entries),1),'color',vcolors{variable},'linewidth',1);
        %plot(xvals,uniqueVarMouseTracked{variable}(:,entries),'color',vcolors{variable}*1.1);
        
    end
    
end

xticks(1:3);
xticklabels({'early','middle','late'});
xtickangle(45);
xlim([0 xvals(end)+1]);
ylim(aloLim);
if variable==1 ylabel('unique R^2'); end
title(vnames(variable));
plot([0.5 xvals(end)+0.5],[0 0],':','color','k');

data=cat(1,uniqueVarMouseTracked{:});
ml=repmat([1;2;3;4;5;6;7;8],3,9);
sl=repmat([1:9],24,1);
dl=repmat([1 1 1 2 2 2 3 3 3],24,1);
gl=cat(1,ones(8,9),2*ones(8,9),3*ones(8,9));
% [p,tbl,stats]=anovan(data(:),{dl(:),gl(:)},'varnames',{'section','variable'},'model','interaction');
% c=multcompare(stats,'dimension',[1 2],'ctype','bonferroni');

end



%% GLM task variables in spatial location

improvCutoff=0.02;
overallCutoff=0.02;
map=redblue(256);
sn='o1d3';

improvement=-(varExp.(sn)-varExp.(sn)(:,1));



%get proportions
categoryMap=zeros(length(improvement),1);
predictive = varExp.(sn)(:,1)>overallCutoff; %more than 5% of variance explained
behavior = improvement(:,3)>improvCutoff; %num licks or lick rate
cueK = improvement(:,2)>improvCutoff; %cue
categoryMap(predictive & behavior,1)=1;
categoryMap(predictive & cueK,1)=2;
categoryMap(predictive & behavior & cueK,1)=3;

catColors=[0.7 0.7 0.7;0.9 0.1 0;0 0.3 0.9;0.99 0 0.99];

mice=unique(roiMouse.(sn));
figure;
for mouse=1:length(mice)
    mn=mice{mouse};
    mouseSel=strcmp(roiMouse.(sn),mn);
    categoryMapMouse=categoryMap(mouseSel);
    load(fullfile(direc,mn,sn,'Fall.mat')); 
    % remove doubles and ROIs on edge of FOV
    iscell = processROIs(iscell,stat,F);
           allROIs = {ones(512,512),ones(512,512),ones(512,512)};
    numROIs=sum(iscell(:,1));
    idx=find(iscell(:,1));
    ROIpixels={};
    for roi=1:numROIs
        roiNum=idx(roi);
        ROIpixels{roi,1}=sub2ind([512 512],stat{roiNum}.ypix,stat{roiNum}.xpix);
        for rgb=1:3
            allROIs{rgb}(ROIpixels{roi,1})=catColors(categoryMapMouse(roi)+1,rgb);
        end
    end

   subplot(3,3,mouse);



    ROIMasks=cat(3,allROIs{:});
    A = ones(512,512);
    A(ROIMasks(:,:,1) == 1) = 0;
    im=imshow(ROIMasks);
    im.AlphaData=A;

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
%% nanste function
%standard error, omitting NaNs
function  ste=nanste(dat,dimension)
ste=nanstd(dat,[],dimension)./sqrt(sum(~isnan(dat),dimension));
end
