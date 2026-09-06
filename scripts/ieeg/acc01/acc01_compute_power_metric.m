function out = acc01_compute_power_metric(pair, events, contract, varargin)
% Optional Name=Value diagnostic mode exposes the already-fit aperiodic
% coefficients without changing the frozen three-input numerical path.
opts=struct("ReturnAperiodicDiagnostics",false);
if ~isempty(varargin)
    assert(mod(numel(varargin),2)==0,"Optional arguments must be Name=Value pairs.");
    for k=1:2:numel(varargin)
        assert(string(varargin{k})=="ReturnAperiodicDiagnostics","Unknown option.");
        opts.ReturnAperiodicDiagnostics=logical(varargin{k+1});
    end
end
warningId='stats:statrobustfit:IterationLimit';warningState=warning('query',warningId);warning('error',warningId);warningCleanup=onCleanup(@()warning(warningState));nonconverged=0;fitCount=0;
fs=contract.sampling.analysis_hz; n=round(contract.power.hann_epoch_s*fs); h=hann(n); f=(0:n-1)'*fs/n; fitix=ismembertol(f,contract.power.aperiodic_fit_bins_hz,1e-10); th=ismembertol(f,contract.power.theta_periodogram_bins_hz,1e-10); y=nan(height(events),1);
if ~opts.ReturnAperiodicDiagnostics
    % Do not refactor this branch: P0-2 callers require bit-identical values.
    for e=1:height(events), c=round(events.label_time_s(e)*fs)+1; a=pair.signal(c+round(contract.windows_s.power_pre(1)*fs):c+round(contract.windows_s.power_pre(1)*fs)+n-1); b=pair.signal(c+round(contract.windows_s.power_post(1)*fs):c+round(contract.windows_s.power_post(1)*fs)+n-1); a=(a-mean(a)).*h; b=(b-mean(b)).*h; y(e)=res(b)-res(a); end
    X=[ones(height(events),1),double(events.z_e),double(events{:,endsWith(events.Properties.VariableNames,"_delta_z")|endsWith(events.Properties.VariableNames,"movie_time_linear_z")|endsWith(events.Properties.VariableNames,"movie_time_quadratic_z")})]; keep=all(isfinite([X y]),2); q=X(keep,:)\y(keep); out=struct('slope',q(2),'dP',y,'n',nnz(keep),'robustfit_iteration_limit_count',nonconverged,'robustfit_fit_count',fitCount);
else
    nEvent=height(events); offsetPre=nan(nEvent,1);offsetPost=offsetPre;exponentPre=offsetPre;exponentPost=offsetPre;rawPre=offsetPre;rawPost=offsetPre;residualDelta=offsetPre;
    for e=1:nEvent
        c=round(events.label_time_s(e)*fs)+1; a=pair.signal(c+round(contract.windows_s.power_pre(1)*fs):c+round(contract.windows_s.power_pre(1)*fs)+n-1); b=pair.signal(c+round(contract.windows_s.power_post(1)*fs):c+round(contract.windows_s.power_post(1)*fs)+n-1);
        a=(a-mean(a)).*h; b=(b-mean(b)).*h;
        [ra,offsetPre(e),exponentPre(e),rawPre(e)]=resWithAperiodic(a);
        [rb,offsetPost(e),exponentPost(e),rawPost(e)]=resWithAperiodic(b);
        y(e)=rb-ra; residualDelta(e)=y(e);
    end
    X=[ones(nEvent,1),double(events.z_e),double(events{:,endsWith(events.Properties.VariableNames,"_delta_z")|endsWith(events.Properties.VariableNames,"movie_time_linear_z")|endsWith(events.Properties.VariableNames,"movie_time_quadratic_z")})]; keep=all(isfinite([X y]),2); q=X(keep,:)\y(keep);
    out=struct('slope',q(2),'dP',y,'n',nnz(keep),'robustfit_iteration_limit_count',nonconverged,'robustfit_fit_count',fitCount, ...
        'aperiodic',struct('offset_pre_log10',offsetPre,'offset_post_log10',offsetPost,'delta_offset_log10',offsetPost-offsetPre, ...
        'exponent_pre',exponentPre,'exponent_post',exponentPost,'delta_exponent',exponentPost-exponentPre, ...
        'raw_theta_pre_log10',rawPre,'raw_theta_post_log10',rawPost,'raw_theta_delta_log10',rawPost-rawPre, ...
        'residual_theta_delta_log10',residualDelta));
end
 function r=res(x)
     p=log10(abs(fft(x)).^2/sum(h.^2));fitCount=fitCount+1;
     try
         q=robustfit(log10(f(fitix)),p(fitix));
     catch ME
         if string(ME.identifier)~=string(warningId),rethrow(ME);end
         nonconverged=nonconverged+1;warning('off',warningId);
         q=robustfit(log10(f(fitix)),p(fitix));warning('error',warningId);
     end
      r=mean(p(th)-(q(1)+q(2)*log10(f(th))));
  end
 function [r,offset,exponent,rawTheta]=resWithAperiodic(x)
      p=log10(abs(fft(x)).^2/sum(h.^2));fitCount=fitCount+1;
      try
          q=robustfit(log10(f(fitix)),p(fitix));
      catch ME
          if string(ME.identifier)~=string(warningId),rethrow(ME);end
          nonconverged=nonconverged+1;warning('off',warningId);
          q=robustfit(log10(f(fitix)),p(fitix));warning('error',warningId);
      end
      offset=q(1);exponent=-q(2);rawTheta=mean(p(th));
      r=mean(p(th)-(q(1)+q(2)*log10(f(th))));
 end
clear warningCleanup
end
