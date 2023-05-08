%% This is the supporting code for paper:
% Zhou, Y, You, J, Kumar, U, Weiss, SA, Bragin, A, Engel, J, et al.
% An approach for reliably identifying high-frequency oscillations and
% reducing false-positive detections. Epilepsia Open. 2022; 7: 674– 686.

% https://onlinelibrary.wiley.com/doi/full/10.1002/epi4.12647

% by Yufeng Zhou 2023

%% I. Preparation
save_path=cd;
sel_path=cd;
inputfiles='test.rhfe';

ch_inds=1:16; % channel number
fs=3000; % sampling frequency
Label=cell(length(ch_inds),1);
for i=1:length(Label)
    label=ch_inds(i);
    Label{i,1}=strcat('ch',num2str(label));
end

nyquist=fs/2;
min_freq=80;
max_freq=520;
filtbound=[min_freq,max_freq];
trans_width=0.01;
filt_order=101;
ffrequencies =[0 (1-trans_width)*filtbound(1) filtbound (1+trans_width)*filtbound(2) nyquist]/nyquist;
idealresponse=[0 0 1 1 0 0];
filterweights=firls(filt_order,ffrequencies,idealresponse);

len=0.25;
offset=0.05;
mode="Normalized";
smoothed="No";
weighted="No";
num_freq=max_freq-min_freq+1;
frex=linspace(min_freq,max_freq,num_freq);


%% II. The main function

rhfe_file=inputfiles;
rhfe_id=fullfile(sel_path,rhfe_file);
f=load(rhfe_id,'-mat');

idoutput=strcat(rhfe_file(1,1:end-5),'_processed.rhfe');
idoutput_save=strcat(fullfile(save_path,idoutput));
s=struct();
s=setfield(s,'st_FileData','str_FileType','edf');
s=setfield(s,'st_FileData','s_Start',[]);
duration=getfield(f,'st_FileData','s_Time');
s=setfield(s,'st_FileData','s_Time',duration);
samples=getfield(f,'st_FileData','s_Samples');
s=setfield(s,'st_FileData','s_Samples',samples);
samplerate=getfield(f,'st_FileData','v_SampleRate');
s=setfield(s,'st_FileData','v_SampleRate',samplerate);
s=setfield(s,'st_FileData','s_NumbRec',16);
s=setfield(s,'st_FileData','v_Labels',Label);
s=setfield(s,'st_FileData','s_Scale',1);
EDF=getfield(f,'st_FileData','st_Custom');
s=setfield(s,'st_FileData','st_Custom',EDF);
unit=getfield(f,'st_FileData','v_AmpScaleRec');
s=setfield(s,'st_FileData','v_AmpScaleRec',unit);
MinMaxRec=getfield(f,'st_FileData','v_MinMaxRec');
s=setfield(s,'st_FileData','v_MinMaxRec',MinMaxRec);
s=setfield(s,'st_FileData','str_SigPath',save_path);
s=setfield(s,'st_FileData','str_FileName',idoutput);
s=setfield(s,'st_FileData','str_SigExt','EDF');
s=setfield(s,'st_FileData','s_error',0);
s=setfield(s,'st_FileData','s_Check',1);
save(idoutput_save,'-struct','s');

for ch_ind=1:length(Label)
    fprintf(['Begin to process channel ',num2str(ch_ind),' of file ','\n']);
    ch_name=Label{ch_ind,1};
    
    s=struct();
    s=setfield(s,ch_name,'st_HFOSetting','s_FreqIni',min_freq);
    s=setfield(s,ch_name,'st_HFOSetting','s_FreqEnd',max_freq);
    s=setfield(s,ch_name,'st_HFOSetting','s_EpochTime',600);
    s=setfield(s,ch_name,'st_HFOSetting','s_RMSWindow',3);
    s=setfield(s,ch_name,'st_HFOSetting','s_RMSThres',5);
    s=setfield(s,ch_name,'st_HFOSetting','s_MinWind',6);
    s=setfield(s,ch_name,'st_HFOSetting','s_MinTime',10);
    s=setfield(s,ch_name,'st_HFOSetting','s_NumOscMin',6);
    s=setfield(s,ch_name,'st_HFOSetting','s_BPThresh',3);
    
    file=getfield(f,ch_name);
    data=getfield(file,'v_Intervals');
    fs=getfield(file,'st_HFOInfo','s_Sampling');
    event_time_ind=getfield(file,'st_HFOInfo','m_IntervLims');
    event_lims=getfield(file,'st_HFOInfo','m_EvtLims');
    event_lims_rel=getfield(file,'st_HFOInfo','m_Rel2IntLims');
    
    
    [Idx,Pks,valid_epoch_peak]=f_peakcheck(filterweights,data,len,fs);
    TF=f_gaussain(data,len,fs,frex,offset,mode,smoothed);
    [Power_cell,M_cell,valid_epoch_pow]=f_noise(TF,len,offset,fs,frex,smoothed);
    clear TF;
    [new_M_cell,Contour_delete,valid_epoch_pow]=f_pow_thresholding(Power_cell,M_cell,valid_epoch_pow);
    clear M_cell;
    [CLC,OLC,valid_epoch_pow]=f_split(new_M_cell,valid_epoch_pow);
    clear new_M_cell
    clear OLC;
    [new_CLC,valid_epoch_lv]=f_regroup(CLC);
    [new_CLC,valid_epoch_lv]=f_valley(new_CLC,valid_epoch_lv);
    [the_Group,B,C]=f_find(new_CLC);
    [valid_epoch_off]=f_off_centered(B);
    valid_overall=cell2mat(valid_epoch_peak).*cell2mat(valid_epoch_pow).*cell2mat(valid_epoch_lv).*cell2mat(valid_epoch_off);
    [C_info,Type]=f_features(B,C,Power_cell,valid_overall,len,offset,fs,frex);
    s=setfield(s,ch_name,'st_HFOInfo','s_Sampling',fs);
    str_ChLabel{1,1}=ch_name;
    s=setfield(s,ch_name,'st_HFOInfo','str_ChLabel',str_ChLabel);
    s=setfield(s,ch_name,'st_HFOInfo','str_DetMethod','Topographical');
    s=setfield(s,ch_name,'st_HFOInfo','s_Chldx',ch_ind);
    
    if ~isempty(valid_overall)
        ind_logic=logical(valid_overall);
        event_time_ind_valid=event_time_ind(ind_logic,1:end);
        s=setfield(s,ch_name,'st_HFOInfo','m_IntervLims',event_time_ind_valid);
        event_lims_valid=event_lims(ind_logic,1:end);
        s=setfield(s,ch_name,'st_HFOInfo','m_EvtLims',event_lims_valid);
        event_lims_rel_valid=event_lims_rel(ind_logic,1:end);
        s=setfield(s,ch_name,'st_HFOInfo','m_Rel2IntLims',event_lims_rel_valid);
        Type_valid=Type(ind_logic);
        s=setfield(s,ch_name,'st_HFOInfo','v_EvType',Type_valid);
        data_valid=data(ind_logic,1);
        s=setfield(s,ch_name,'v_Intervals',data_valid);
        s=setfield(s,ch_name,'delete_info','layer_I',valid_epoch_peak);
        s=setfield(s,ch_name,'delete_info','layer_II',valid_epoch_pow);
        s=setfield(s,ch_name,'delete_info','layer_III',valid_epoch_lv);
        s=setfield(s,ch_name,'delete_info','layer_IV',valid_epoch_off);
        s=setfield(s,ch_name,'delete_info','final_result',valid_overall);
    else
        s=setfield(s,ch_name,'st_HFOInfo','m_IntervLims',[]);
        s=setfield(s,ch_name,'st_HFOInfo','m_EvtLims',[]);
        s=setfield(s,ch_name,'st_HFOInfo','m_Rel2IntLims',[]);
        s=setfield(s,ch_name,'st_HFOInfo','v_EvType',[]);
        s=setfield(s,ch_name,'v_Intervals',{});
        s=setfield(s,ch_name,'delete_info','layer_I',valid_epoch_peak);
        s=setfield(s,ch_name,'delete_info','layer_II',valid_epoch_pow);
        s=setfield(s,ch_name,'delete_info','layer_III',valid_epoch_lv);
        s=setfield(s,ch_name,'delete_info','layer_IV',valid_epoch_off);
        s=setfield(s,ch_name,'delete_info','final_result',valid_overall);
    end
    save(idoutput_save,'-struct','s','-append');
end


function [Idx,Pks,valid_epoch_peak]=f_peakcheck(filterweights,data,len,fs)
valid_epoch_peak=cell(length(data),1);
Idx=cell(length(data),1);
Pks=cell(length(data),1);
for i=1:length(data)
    d=transpose(data{i,1});
    mid_point=floor(length(d)/2)+1;
    d=d(1,mid_point-len*fs:mid_point+len*fs);
    d_filt=filtfilt(filterweights,1,d);
    peakThresh1=mean(d_filt)+2*std(d_filt);
    peakThresh2=mean(d_filt)-2*std(d_filt);
    locs1=find((d_filt>=peakThresh1)==1);
    peaks1=d_filt(locs1);
    locs2=find((d_filt<=peakThresh2)==1);
    peaks2=d_filt(locs2);
    idx=[locs1,locs2];
    pks=[peaks1,peaks2];
    start_ind=fs*(len-0.1);
    l1=(idx>=start_ind);
    end_ind=fs*(2*len-0.1);
    l2=(idx<=end_ind);
    res=sum(l1.*l2);
    if res/length(idx)>0.60
        valid_epoch_peak{i,1}=1
    else
        valid_epoch_peak{i,1}=0;
    end
    Idx{i,1}=idx;
    Pks{i,1}=pks;
end
end

function TF=f_gaussain(data,len,fs,frex,offset,mode,smoothed)
if isempty(data)
    TF={};
else
    TF=cell(size(data,1),1);
    for i=1:size(data,1)
        d=cell2mat(data(i,:))./1000;
        middle_point=floor(size(d,1)/2)+1;
        d=d(middle_point-len*fs:middle_point+len*fs);
        dX=fft(d,length(d));
        cycles=3;
        tvec=-len:1/fs:len;
        half_len=floor(length(d)/2)+1;
        wvec=(2*pi/length(d))*(0:(length(d)-1))*fs;
        wvec_half=wvec(1:half_len);
        tf=zeros(length(frex),length(tvec));
        for j=1:length(frex)
            sigma=(1/frex(j))*cycles;
            gau_win_X=zeros(length(d),1);
            gau_win_X(1:half_len,1)=exp(-0.5.*realpow(wvec_half-2*pi*frex(j),2).*(sigma.^2));
            gau_win_X=gau_win_X.*sqrt(length(d))./norm(gau_win_X,2);
            tf(j,1:end)=ifft(dX.*gau_win_X)./sqrt(sigma);
        end
        tf=abs(tf);
        tf=tf(1:end,1+offset*fs:end-offset*fs);
        if mode=="Normalized"
            tf_sum=sum(tf,2);
            for k=1:size(tf,1)
                tf(k,:)=tf(k,:)/tf_sum(k,1);
            end
        elseif mode=="Raw"
            fprintf('Remain unchanged\n');
        end
        if smoothed=="Yes"
            smooth_tf=zeros(size(tf));
            smooth_steps=10;
            for k=1:smooth_steps
                idx=360*ones(smooth_steps,1);
                temp=movmean(tf,idx(k,1),2);
                if k<=5
                    temp=movmean(temp,idx(k,1)./6,1);
                end
                smooth_tf=smooth_tf+temp;
            end
            TF{i,1}=smooth_tf./smooth_steps;
            fprintf(['The process of file number ',num2str(i),' is finished\n']);
        elseif smoothed=="No"
            TF{i,1}=tf;
        end
    end
    fprintf('Spetrum calulation is finised!!!\n');
end
end

function [Power_cell,M_cell,valid_epoch_pow]=f_noise(TF,len,offset,fs,frex,smoothed)
if isempty(TF)
    Power_cell={};
    M_cell={};
    valid_epoch_pow={};
else
    Power_cell=cell(size(TF,1),1);
    M_cell=cell(size(TF,1),1);
    valid_epoch_pow=cell(size(TF,1),1);
    if smoothed=="Yes"
        s1=1;
        s2=0.7;
    elseif smoothed=="No"
        s1=2;
        s2=0.7;
    end
    fprintf('Begin the noise cancellation !!!\n');
    for i=1:size(TF,1)
        pm=TF{i,1};
        epoch_len=0.02*fs;
        num_epochs=(size(pm,2)-1)/epoch_len;
        pow_thresh=mean(pm,'all')+s1*std(pm,1,'all');
        cnt_thresh=s2*length(frex);
        j=1;
        while j<=num_epochs
            pow_val=mean(pm(:,(j-1)*epoch_len+1:j*epoch_len),'all');
            pow_cond=(pow_val>=pow_thresh);
            cnt_val=sum(mean(pm(:,(j-1)*epoch_len+1:j*epoch_len),2)>mean(TF{i,1},2));
            cnt_cond=(cnt_val>=cnt_thresh);
            cond=pow_cond*cnt_cond;
            if cond==1
                valid_epoch_pow{i,1}=0;
                j=100000;
            else
                valid_epoch_pow{i,1}=1;
                j=j+1;
            end
        end
        x=linspace(-(len-offset),(len-offset),2*(len-offset)*fs+1);
        y=transpose(frex);
        [X,Y]=meshgrid(x,y);
        Z=pm;
        M=contour(X,Y,Z,20);
        Power_cell{i,1}=Z;
        M_cell{i,1}=M;
    end
    fprintf('Finish the noise cancellation process!!!\n');
end
end

function [new_M_cell,Contour_delete,valid_epoch_pow]=f_pow_thresholding(Power_cell,M_cell,valid_epoch_pow)
if isempty(Power_cell)
    new_M_cell={};
    Contour_delete={};
    valid_epoch_pow={};
else
    new_M_cell={};
    Contour_delete=zeros(size(M_cell,1),1);
    fprintf('Begin the power thresholding process!!!\n');
    for i=1:size(M_cell,1)
        Z=Power_cell{i,1};
        M=M_cell{i,1};
        power_ind=zeros(10000000,1);
        cnt=1;
        position=1;
        while position<size(M,2)
            power_ind(cnt,1)=position;
            cnt=cnt+1;
            num_vertices=M(2,position);
            position=position+num_vertices+1;
        end
        power_ind(power_ind==0)=[];
        new_M={};
        contour_delete=0;
        flag=0;
        power_threshold=0.2*(max(Z,[],'all')-min(Z,[],'all'))+min(Z,[],'all');
        for j=1:length(power_ind)-1
            power_level=M(1,power_ind(j));
            if power_level<power_threshold
                contour_delete=contour_delete+1;
                continue
            else
                flag=flag+1;
                ind_start=power_ind(j,1);
                ind_end=power_ind(j+1,1)-1;
                new_M{flag}=M(:,ind_start:ind_end);
            end
        end
        new_M_cell{i}=new_M;
        Contour_delete(i,1)=contour_delete;
        if isempty(new_M)
            valid_epoch_pow{i,1}=0*valid_epoch_pow{i,1};
        else
            valid_epoch_pow{i,1}=1*valid_epoch_pow{i,1};
        end
    end
    fprintf('Finish the power thresholding process!!!\n');
end
end

function [CLC,OLC,valid_epoch_pow]=f_split(new_M_cell,valid_epoch_pow)
if isempty(new_M_cell)
    CLC={};
    OLC={};
    valid_epoch_pow={};
else
    CLC={};
    OLC={};
    fprintf('Start the contour classfication!!!\n');
    for i=1:size(new_M_cell,2)
        new_M=new_M_cell{i};
        clc={};
        olc={};
        ind_clc=0;
        ind_olc=0;
        for j=1:size(new_M,2)
            contour=new_M{j};
            if (contour(1,2)==contour(1,end))&&(contour(2,2)==contour(2,end))
                ind_clc=ind_clc+1;
                clc{ind_clc}=contour;
            else
                ind_olc=ind_olc+1;
                olc{ind_olc}=contour;
            end
        end
        CLC{i}=clc;
        OLC{i}=olc;
        if isempty(clc)
            valid_epoch_pow{i,1}=0*valid_epoch_pow{i,1};
        else
            valid_epoch_pow{i,1}=1*valid_epoch_pow{i,1};
        end
    end
    fprintf('Finish the classificaton process!!!\n');
end
end

function [new_CLC,valid_epoch_lv]=f_regroup(CLC)
if isempty(CLC)
    new_CLC={};
    valid_epoch_lv={};
else
    valid_epoch_lv=cell(size(CLC,1),1);
    new_CLC={};
    parfor i=1:size(CLC,2)
        unassign_clcs=CLC{i};
        new_clc={};
        cnt=1;
        while length(unassign_clcs)>3
            max_area=-10000;
            max_area_ind=0;
            for j=1:length(unassign_clcs)
                clc=unassign_clcs{j};
                area=polyarea(clc(1,2:end),clc(2,2:end));
                if area>max_area
                    max_area=area;
                    max_area_ind=j;
                end
            end
            b=unassign_clcs{max_area_ind};
            group_ind=[];
            group_ind(length(group_ind)+1)=max_area_ind;
            for j=1:length(unassign_clcs)
                clc=unassign_clcs{j};
                b_left=min(b(1,2:end));
                x_left=min(clc(1,2:end));
                cond1=1*(b_left<x_left);
                b_right=max(b(1,2:end));
                x_right=max(clc(1,2:end));
                cond2=1*(b_right>x_right);
                b_bottom=min(b(2,2:end));
                x_bottom=min(clc(2,2:end));
                cond3=1*(b_bottom<x_bottom);
                b_top=max(b(2,2:end));
                x_top=max(clc(2,2:end));
                cond4=1*(b_top>x_top);
                poly1=polyshape(b(1,2:end),b(2,2:end));
                poly2=polyshape(clc(1,2:end),clc(2,2:end));
                inside=isinterior(poly1,poly2.Vertices);
                cond5=(sum(inside)==size(inside,1));
                cond=cond1*cond2*cond3*cond4*cond5;
                if cond==1
                    group_ind(length(group_ind)+1)=j
                end
            end
            if length(group_ind)>3
                new_clc{cnt}=unassign_clcs(1,group_ind);
                cnt=cnt+1;
            end
            unassign_clcs(:,group_ind)=[];
        end
        new_CLC{i}=new_clc;
        fprintf(['The regrouping of number ',num2str(i),' epoch is finished\n']);
        if isempty(new_clc)
            valid_epoch_lv{i,1}=0;
        else
            valid_epoch_lv{i,1}=1;
        end
    end
end
end

function [new_CLC,valid_epoch_lv]=f_valley(new_CLC,valid_epoch_lv)
if isempty(new_CLC)
    new_CLC={};
    valid_epoch_lv={};
else
    valley_clc=cell(1,size(new_CLC,2));
    for i=1:size(new_CLC,2)
        cnt=0;
        clc_epoch=new_CLC{i};
        for j=1:length(clc_epoch)
            clc_group=clc_epoch{j};
            power_max=-10000;
            power_max_ind=0;
            for k=1:length(clc_group)
                clc=clc_group{k};
                power=clc(1,1);
                if power>power_max
                    power_max=power;
                    power_max_ind=k;
                end
            end
            c=clc_group{power_max_ind};
            for k=1:length(clc_group)
                clc=clc_group{k};
                c_left=min(c(1,2:end));
                x_left=min(clc(1,2:end));
                cond1=1*(c_left<x_left);
                c_right=max(c(1,2:end));
                x_right=max(clc(1,2:end));
                cond2=1*(c_right>x_right);
                c_bottom=min(c(2,2:end));
                x_bottom=min(clc(2,2:end));
                cond3=1*(c_bottom<x_bottom);
                c_top=max(c(2,2:end));
                x_top=max(clc(2,2:end));
                cond4=1*(c_top>x_top);
                cond=cond1*cond2*cond3*cond4
                if cond==1
                    cnt=cnt+1;
                    coords=[i,j];
                    valley_clc{1,i}{1,cnt}=coords
                end
            end
        end
        fprintf(['The valley checking of ',num2str(i),' epoch is finished\n']);
    end
    for i=1:length(valley_clc)
        valley_clc_group_info=valley_clc{i};
        if isempty(valley_clc_group_info)==0
            for j=1:length(valley_clc_group_info)
                info=valley_clc_group_info{j};
                event_ind=info(1); % the index of event
                group_ind=info(2); % the index of group in the event
                new_CLC{1,event_ind}{1,group_ind}=[]; % delete valleys (clcs in event)
            end
        end
    end
    for i=1:length(new_CLC)
        clcs=new_CLC{i}(~cellfun('isempty',new_CLC{i}));
        new_CLC{i}=clcs;
        if isempty(clcs)
            valid_epoch_lv{i,1}=0*valid_epoch_lv{i,1};
        else
            valid_epoch_lv{i,1}=1*valid_epoch_lv{i,1};
        end
    end
end
end

function [the_Group,B,C]=f_find(new_CLC)
if isempty(new_CLC)
    the_Group={};
    B={};
    C={};
else
    the_Group=cell(length(new_CLC),1);
    B=cell(length(new_CLC),1);
    C=cell(length(new_CLC),1);
    fprintf('Strat to extract the Boundary and Peak!!!\n');
    for i=1:length(new_CLC)
        clc_epoch=new_CLC{i};
        if isempty(clc_epoch)
            the_Group{i,1}={};
            B{i,1}={};
            C{i,1}={};
        else
            max_power_level=zeros(length(clc_epoch),1);
            max_power_ind=zeros(length(clc_epoch),1);
            Cs=cell(length(clc_epoch),1);
            Bs=cell(length(clc_epoch),1);
            for j=1:length(clc_epoch)
                clc_group=clc_epoch{j};
                size_max=-1000000;
                size_max_ind=0;
                power_max=-1000000;
                power_max_ind=0;
                power_max_size=0;
                for l=1:length(clc_group)
                    clc=clc_group{l};
                    clc_size=polyarea(clc(1,2:end),clc(2,2:end));
                    clc_power=clc(1,1);
                    if clc_size>size_max
                        size_max=clc_size;
                        size_max_ind=l;
                    end
                    if clc_power>power_max
                        power_max=clc_power;
                        power_max_ind=l;
                        power_max_size=clc_size;
                    elseif clc_power==power_max
                        if clc_size>power_max_size
                            power_max=clc_power;
                            power_max_ind=l;
                            power_max_size=clc_size;
                        end
                    end
                end
                Bs{j,1}=clc_group{size_max_ind};
                max_power_level(j,1)=power_max;
                max_power_ind(j,1)=power_max_ind;
                Cs{j,1}=clc_group{power_max_ind};
            end
            [~,max_group]=max(max_power_level);
            the_Group{i,1}=clc_epoch{max_group};
            C{i,1}=Cs{max_group,1};
            B{i,1}=Bs{max_group,1};
        end
    end
    fprintf('Finish the extracting process!!!\n');
end
end

function [valid_epoch_off]=f_off_centered(B)
if isempty(B)
    valid_epoch_off={};
else
    valid_epoch_off=cell(length(B),1);
    fprintf('Start the off-center checking!!!\n');
    for i=1:length(B)
        b=B{i,1};
        if isempty(b)==0
            xmin=min(b(1,2:end));
            xmax=max(b(1,2:end));
            if (xmax>=-0.05)*(xmin<=0.05)==1
                valid_epoch_off{i,1}=1;
            else
                valid_epoch_off{i,1}=0;
            end
        else
            valid_epoch_off{i,1}=0;
        end
    end
    fprintf('Finish the off-center checking process!!!\n');
end
end

function [C_info,Type]=f_features(B,C,Power_cell,valid_overall,len,offset,fs,frex)
if isempty(B)
    C_info=[];
    Type=[];
else
    Power_mean=zeros(size(B,1),1);
    Power_sum=zeros(size(B,1),1);
    Freq_weighted=zeros(size(B,1),1);
    Freq_sum=zeros(size(B,1),1);
    Type=zeros(size(B,1),1);
    for i=1:size(C,1)
        if valid_overall(i,1)==1
            Z=Power_cell{i,1};
            c=C{i,1};
            xv=c(1,2:end);
            yv=c(2,2:end);
            x=linspace(-(len-offset),(len-offset),2*(len-offset)*fs+1);
            y=transpose(frex);
            [X,]=meshgrid(x,y);
            xq=transpose(reshape(X,[size(X,1)*size(X,2),1]));
            yq=repmat(frex,1,length(x));
            [in,on]=inpolygon(xq,yq,xv,yv);
            grid_ind=zeros(1,length(xq(in))+length(xq(on)));
            grid_ind(1:length(xq(in)))=find(in);
            grid_ind(length(xq(in))+1:end)=find(on);
            power_sum=0;
            freq_sum=0;
            for j=1:size(grid_ind,2)
                ind=grid_ind(j);
                z_x=xq(ind);
                z_y=yq(ind);
                power=Z(frex==z_y,x==z_x);
                freq_sum=freq_sum+z_y*power;
                power_sum=power_sum+power;
            end
            power_mean=power_sum/size(grid_ind,2);
            Power_mean(i,1)=power_mean;
            Power_sum(i,1)=power_sum;
            freq_weighted=freq_sum/power_sum;
            Freq_weighted(i,1)=freq_weighted;
            Freq_sum(i,1)=freq_sum;
            if freq_weighted<=240
                Type(i,1)=1;
            elseif freq_weighted>240
                Type(i,1)=2;
            end
            fprintf(['Quantifying in event ',num2str(i),' is finished\n']);
        end
    end
    C_info=[Power_mean,Power_sum,Freq_weighted,Freq_sum];
end
end
