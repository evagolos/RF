function DISP = ad_SACPZ(raw,zeros,poles,gain,samprate)
% remove instrument response from raw data into DISP
% the poles, zeros and gain values should get from SAC pole zero files (DISP)
% figure(101) shows RESP ampitude in VEL and DISP 
% pylin.patty 2014.03

% should write one for RESP file, and normalized RESP amplitude at sensitive frequency and then multiply gain.
% pylin.patty 2014.04

%isfigure = 0;
%raw=flat_hanning_win(1:length(raw),raw,1,length(raw),200);

N=length(raw);
T=N*(1./samprate);

% way 1 
if mod(N,2)
     faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
else
     faxis = [0:N/2,-N/2+1:-1]*(1/T);
end
w = faxis.*2*pi;
resp = ones(size(w));
for ip = 1:length(poles)
        resp = resp./(1i*w - poles(ip));
end
for ip = 1:length(zeros)
        resp = resp.*(1i*w - zeros(ip));
end

%resp=resp;
%index= knnsearch(faxis',0.2);
resp = resp*gain;

% plotting ===
%figure(101)
%clf
%subplot(2,1,2)
%loglog(faxis,abs(resp));hold on;grid on;
%title('DISP response');
%subplot(2,1,1)
%loglog(faxis,abs(resp)./faxis/2/pi);grid on;
%title('VEL response');
% =============
npoles = length(poles);

lo_corner = 0.009;
lo_w=2*pi*lo_corner;
hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );

% hi_corner = 20;
% hi_w=2*pi*hi_corner;
% hpfiltfrq=  hpfiltfrq - ( ((w./hi_w).^(2*npoles))./(1+(w./hi_w).^(2*npoles)) );


norm_trans=hpfiltfrq./resp;    % this is normalization transfer function
norm_trans(find(isnan(norm_trans))) = 0;


% if ~isfigure
%     figure(101)
%     clf
%     %plot(w/2.0/pi,hpfiltfrq,'bx');
%     plot(w/2.0/pi,abs(norm_trans),'bx');
% %pause
% %    %xlim([-10 10]);
% end
%close(101)



fftdata = raw;
fftdata = fftdata(:)./norm_trans(:);
fftdata(isnan(fftdata)==1)=0;
DISP = real(ifft(fftdata));
DISP = detrend(DISP);
DISP = DISP-mean(DISP);







return
