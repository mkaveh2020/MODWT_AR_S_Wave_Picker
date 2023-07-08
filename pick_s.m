clc
clear
close all

%%
wname ='db2';
%%                                                                                                                                                                                                                                                                                                                      


[FileName,PathName] = uigetfile('*.SAC','Open ...','MultiSelect', 'on');

File = strcat(PathName,FileName);

%][pytre
 %  iytr
for i = 1 :max(size(File))
    s=readsac(File{i});
    tttt=((0:(length(s.DATA1)-1))*s.DELTA)+s.B;
    
    if (s.CMPAZ==0 && s.CMPINC==0)
        %plot(t,s.DATA1)
        tr33=s.DATA1;
        noise = randn(size(tr33))*rms(tr33)
        tr33 = tr33 + .3* noise
        %title('BHZ')
    else if (s.CMPAZ==0 && s.CMPINC==90)
            tr22=s.DATA1;
            noise = randn(size(tr22))*rms(tr22)
            tr22 = tr22 + .3* noise
            % plot(t,s.DATA1)
            %title('BHN')
        else if (s.CMPAZ==90 && s.CMPINC==90)
                tr11=s.DATA1;
                noise = randn(size(tr11))*rms(tr11)
                tr11 = tr11 + .3*noise
                %plot(t,s.DATA1)
                %title('BHE')
            end
        end
    end
    %axis([300 600 -2000 2000])
    %  zoom on;
    %  %pause();
    %  zoom off;
    % [x,~]=ginput(1);
    % x1 = [x,x];
    % y1 = [max(abs(s.DATA1)),-max(abs(s.DATA1))];
    % hold on
    % plot(x1,y1,'r')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs=100;
% lpFilt = designfilt('lowpassfir','PassbandFrequency',0.05/Fs, ...
%          'StopbandFrequency',5/Fs,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',65,'DesignMethod','kaiserwin');
    lpFilt=  designfilt('lowpassfir', 'PassbandFrequency', 0.25, ...
                 'StopbandFrequency',2.5, ...
                 'PassbandRipple', 1, 'StopbandAttenuation', 70, ...
                 'SampleRate', Fs,'DesignMethod','kaiserwin');
             
rp = 5;           % Passband ripple
rs = 70;          % Stopband ripple
fs = 100;        % Sampling frequency
f = [2.5 3.5];    % Cutoff frequencies
a = [1 0];        % Desired amplitudes
dev = [(10^(rp/20)-1)/(10^(rp/20)+1)  10^(-rs/20)]; 
[n,fo,ao,w] = firpmord(f,a,dev,fs);
b = firpm(n,fo,ao,w);
freqz(b,1,100,fs)            
  
tr11=tr11-mean(tr11);
tr1=filter(b,1,tr11);
tr22=tr22-mean(tr22);
tr2=filter(b,1,tr22);
tr33=tr33-mean(tr33);


tr3=filter(b,1,tr33);
   
% tr1=filtfilt(lpFilt,tr11);
% tr2=filtfilt(lpFilt,tr22);
% tr3=filtfilt(lpFilt,tr33);

tr1 = tr1 - mean(tr1);
tr1 = detrend(tr1,'linear');
tr1 = tukeywin(length(tr1),0.01).*tr1;
tr2 = tr2 - mean(tr2);
tr2 = detrend(tr2,'linear');
tr2 = tukeywin(length(tr2),0.01).*tr2;
tr3 = tr3 - mean(tr3);
tr3 = detrend(tr3,'linear');
tr3 = tukeywin(length(tr3),0.01).*tr3;% taparing the signal
%remove noise for real data
%tr = tr + 0.0*rms(tr)*randn(size(tr));
%  figure
%  plot(t,tr)
%% ba data removal of first and end of data
tr1=tr1(100:length(tr1)-100);
tr2=tr2(100:length(tr2)-100);
tr3=tr3(100:length(tr3)-100);

k=1;   %number of your level wavelet
j=1;   %only 1 or 2 you can choise
W_tr1 = modwt(tr1',wname);
W_tr2 = modwt(tr2',wname);
W_tr3 = modwt(tr3',wname);
% W_tr1=(W_tr1+W_tr1.*sign(W_tr1))/2;
% W_tr2=(W_tr2+W_tr2.*sign(W_tr2))/2;
% W_tr3=(W_tr3+W_tr3.*sign(W_tr3))/2;

%%
%read the seismogram


% [t,tr] = readsac('TU.KSN.U.SAC');
% tr = detrend(tr,'linear');
% tr = tr - mean(tr);
% tr = tr + 0.0*rms(tr)*randn(size(tr));
% figure
% plot(t,tr)
%%
%wavelet transform using wname==db4
% W_tr = modwt(tr',wname);

%%
%envelope estimation using hilbert transform
env1 = zeros(size(W_tr1,1)-1,size(W_tr1,2));


for i=1:size(W_tr1,1)-1
    env1(i,:) = envelope(W_tr1(i,:),50,'peak');
    env2(i,:) = envelope(W_tr2(i,:),50,'peak');
    env3(i,:) = envelope(W_tr3(i,:),50,'peak');    
end



% peak_env1=zeros(size(env1));
% for i=1:size(env1,1)
%     [peaks,locs]= findpeaks(env1(i,:),'MinPeakHeigh',.9*max(env1(i,:)));
%     peak_env1(i,locs)=ones(size(peaks));
% end
% sum_env1=sum(peak_env1);
% 
% subplot(221);plot(sum_env1)
% hold on
% subplot(222);plot(tr2)
% subplot(2,2,3:4);imagesc(env1)
LL=length(tr11);
[envmax1 btest1]=find(env1==max(max(env1(3:9,fix(LL/3):fix(2*LL/3)))));
[envmax2 btest2]=find(env2==max(max(env2(3:9,fix(LL/3):fix(2*LL/3)))));
[envmax3 btest3]=find(env3==max(max(env3(3:9,fix(LL/3):fix(2*LL/3)))));
qq1=(env1(envmax1,:));
[a b1]=find(qq1>=0.5*max(qq1));
b1=b1(1);
[r rr]=zcdd(diff(qq1));
r=r.*(r<b1-1);
[t y]=min(abs(r-b1));
b1=r(y);
qq2=(env2(envmax2,:));
[a b2]=find(qq2>=0.5*max(qq2));
b2=b2(1);
[r rr]=zcdd(diff(qq2));
r=r.*(r<b2-1);
[t y]=min(abs(r-b2));
b2=r(y);
qq3=(env3(envmax3,:));
[a b3]=find(qq3>=0.5*max(qq3));
b3=b3(1);
[r rr]=zcdd(diff(qq3));
r=r.*(r<b3-1);
[t y]=min(abs(r-b3));
b3=r(y);
% b=fix((b1+b2+b3)/3);
% qq1=(env1(5,:));
% [a b1]=find(qq1>0.455*max(env1(:)));
% qq2=(env2(5,:));
% [a b2]=find(qq2>0.455*max(env2(:)));
% qq3=(env3(5,:));
% [a b3]=find(qq3>0.455*max(env3(:)));
% b=fix((b1(1)+b2(1)+b3(1))/3);



figure;plot(tr1)
hold on;plot([b1,b1],[max(abs(tr1)),-max(abs(tr1))],'k--','linewidth',1)
figure;plot(tr2)
hold on;plot([b2,b2],[max(abs(tr2)),-max(abs(tr2))],'r--','linewidth',1)
figure;plot(tr3)
hold on;plot([b3,b3],[max(abs(tr3)),-max(abs(tr3))],'b--','linewidth',1)
% qqvar2(1)=1/(var(tr1(b1:end))/var(tr1(1:b1)));
% qqvar2(2)=1/(var(tr2(b2:end))/var(tr2(1:b2)));
% qqvar2(3)=1/(var(tr3(b3:end))/var(tr3(1:b3)));
% qqvar(1)=qqvar2(1)/sum(qqvar2);
% qqvar(2)=qqvar2(2)/sum(qqvar2);
% qqvar(3)=qqvar2(3)/sum(qqvar2);
% 
% [a4 b4]=min(qqvar);
% figure;plot(tr2)
% hold on;plot([(qqvar(1)*b1+qqvar(2)*b2+qqvar(3)*b3),(qqvar(1)*b1+qqvar(2)*b2+qqvar(3)*b3)],[max(abs(tr2)),-max(abs(tr2))],'r.-','linewidth',2)


% figure
% imagesc(env1)
% set(gca,'YDir','normal')
%
% env2 = zeros(size(W_tr2,1)-1,size(W_tr2,2));
% a=size(W_tr2);
%
% for i=1:size(W_tr2,1)-1
%     env2(i,:) = envelope(W_tr2(i,:));
% end
%  peak_env2=zeros(size(env1));
%     for i=1:size(env2,1)
%         [peaks,locs]= findpeaks(env2(i,:),'MinPeakHeigh',.5*max(env2(i,:)));
%         peak_env2(i,locs)=peaks;
%     end
%
% sum_env2=sum(peak_env2);
% figure
% imagesc(env2)
% set(gca,'YDir','normal')
%
% env3 = zeros(size(W_tr3,1)-1,size(W_tr3,2));
% a=size(W_tr3);
%
% for i=1:size(W_tr3,1)-1
%     env3(i,:) = envelope(W_tr3(i,:));
% end
%  peak_env3=zeros(size(env3));
%     for i=1:size(env3,1)
%         [peaks,locs]= findpeaks(env3(i,:),'MinPeakHeigh',.5*max(env3(i,:)));
%         peak_env3(i,locs)=peaks;
%     end
%     sum_env3=sum(peak_env3);
% figure
% imagesc(env3)
% set(gca,'YDir','normal')

% figure
% imagesc(env)
%%
% %plot the envelop sum
% % figure
% DF1=(abs(env1));
% DF2=(abs(env2));
% DF3=(abs(env3));
% % plot(t,DF)
% %%
% %define the window for samples ER1
% l = 100;
% DF_new1=[zeros(1,l),DF1,zeros(1,l)];
% DF_new2=[zeros(1,l),DF2,zeros(1,l)];
% DF_new3=[zeros(1,l),DF3,zeros(1,l)];
%
% ER1 = zeros(size(DF1));
% ER2 = zeros(size(DF2));
% ER3 = zeros(size(DF3));
% aa=.1:.1:1;
% % mm=length(aa);
% % for j=1:mm
% %     nn(j)=aa(j);
%
% for i = (l+1):(length(DF1)+l)
%     ER1(i-l)=sum(DF_new1(i:i+l))/sum(DF_new1(i-l:i));
% end
% ER1(1:l)=ER1(l+1)*ones(1,l);
% ER1(end-l+1:end)=ER1(end-l-2)*ones(1,l);
%
%
% for i = (l+1):(length(DF2)+l)
%     ER2(i-l)=sum(DF_new2(i:i+l))/sum(DF_new2(i-l:i));
% end
% ER2(1:l)=ER2(l+1)*ones(1,l);
% ER2(end-l+1:end)=ER2(end-l-2)*ones(1,l);
%
%
% for i = (l+1):(length(DF3)+l)
%     ER3(i-l)=sum(DF_new3(i:i+l))/sum(DF_new3(i-l:i));
% end
% ER3(1:l)=ER3(l+1)*ones(1,l);
% ER3(end-l+1:end)=ER3(end-l-2)*ones(1,l);
%
% w_peaks1 = diff(ER1);
% w_peaks1 = w_peaks1.*(w_peaks1>0);
%
% w_peaks2 = diff(ER2);
% w_peaks2 = w_peaks2.*(w_peaks2>0);
%
% w_peaks3 = diff(ER3);
% w_peaks3 = w_peaks3.*(w_peaks3>0);
%
%
% [peak1,tt1]=findpeaks(env1,'MinPeakHeight',.2,'MinPeakDistance',l);%,'Threshold',.1);
% [peak2,tt2]=findpeaks(env2,'MinPeakHeight',.2,'MinPeakDistance',l);%,'Threshold',.1);
% [peak3,tt3]=findpeaks(env3,'MinPeakHeight',.2,'MinPeakDistance',l);%,'Threshold',.1);
% % []
%
%
% %%
% %calculating ER2
% ER2=ER1.*abs(DF);
% % figure;plot(t,ER2,'k')
% figure
% % subplot(3,2,1:2)
% figure;plot(t,tr,'k')
% if (wname(1:3)=='db4')
%     title('TU.KSN.U.SAC ')
% elseif(wname(1:3)=='db2')
%     title('TU.KSN.U.SAC ')
% elseif (wname(1:3) =='db3')
%     title('Dubichease3')
% elseif (wname(1:3)=='db1')
%     title ('Dubichease 1')
% elseif (wname(1:3)=='db8')
%     title ('dubicheas 8')
% elseif (wname=='sym2')
%     title ('symlet 2 "TU.KSN.U.SAC"')
%     elseif (wname =='sym3')
%     title ('symlet 3')
%      elseif (wname =='sym8')
%     title ('symlet 8')
% end
% % zoom on;
% % %pause();
% % zoom off;
% %%
% hold on
%   plot([t(tt),t(tt)],[max(abs(tr)),-max(abs(tr))],':k','linewidth',1)
%
%  tt1= input('tr ra vared konid:');
%  hold on
%  plot([tt1,tt1],[max(abs(tr)),-max(abs(tr))],'k--','linewidth',1)
% % subplot(3,2,3)
% figure;imagesc(env)
% colormap('gray')
% set(gca,'YDir','normal')
% % subplot(3,2,4)
% figure;plot(t,DF,'k')
% % subplot(3,2,5:6)
% figure;plot(w_peaks,'k')
% % hold on
% % plot(diff(ER1),'r')
% title(sprintf('MinPeakHeight=%0.1f',.4))
% % subplot(3,2,6)
% % plot(t,ER2)
% hold on;
% plot(tt,peak,'*k')
% saveas(gcf,sprintf('MinPeakHeight_%0.1f_%s.jpg',.4,wname),'jpeg')
% Fs=100;
% lpFilt = designfilt('lowpassfir','PassbandFrequency',0.05/Fs, ...
%          'StopbandFrequency',5/Fs,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',65,'DesignMethod','kaiserwin');
% aa2=filter(lpFilt,aa);
% fvtool(lpFiilt);

% 
% aa2=tr2;
% aa2=filtfilt(lpFilt,aa2);
% 
% L=length(aa2);
% bb=fft(aa2.*hamming(length(aa2)));
% P2 = abs(bb/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% figure;plot(f,P1) 
% 
bL1=min(b1-1,500);bH1=min(500,length(tr11)-b1);
bL2=min(b2-1,500);bH2=min(500,length(tr11)-b2);
bL3=min(b3-1,500);bH3=min(500,length(tr11)-b3);

qualityfvar(1)=var(tr11(b1-bL1:b1+bH1))+var(tr11(b2-bL2:b2+bH2))+var(tr11(b3-bL3:b3+bH3));
qualityfvar(2)=var(tr22(b1-bL1:b1+bH1))+var(tr22(b2-bL2:b2+bH2))+var(tr22(b3-bL3:b3+bH3));
qualityfvar(3)=var(tr33(b1-bL1:b1+bH1))+var(tr33(b2-bL2:b2+bH2))+var(tr33(b3-bL3:b3+bH3));
% qualityfvar(1)=var(env1(:));
% qualityfvar(2)=var(env2(:));
% qualityfvar(3)=var(env3(:));
[op we]=max(qualityfvar);
if we==1
    fsig=tr11;
    bf=b1;
elseif we==2
    fsig=tr22;
    bf=b2;
else
    fsig=tr33;
    bf=b3;
end



%% ARMA
cu=0;
qqp=-1000:20:0;
for i=-1000:20:0
    cu=cu+1;
    y2=fsig(bf+i-100:bf+i+100);
    y2=y2/(max(y2));
    z = iddata(y2);
    m2 = ar(z,180,'yw');   
[yyh fitee(cu) yyth2]=compare(z,m2);
end
[ww tt]=min(fitee);
figure;plot(fsig)
hold on;plot([bf+qqp(tt),bf+qqp(tt)],[max(abs(fsig)),-max(abs(fsig))],'r-','linewidth',1.5)
pause;
[tt11,~]= ginput(1);
fprintf(1,'difference = %f\n',round(tt11-(bf+qqp(tt)))*s.DELTA);
                  
hold on;plot([tt11,tt11],[max(abs(fsig)),-max(abs(fsig))],'k--','linewidth',1.5)
