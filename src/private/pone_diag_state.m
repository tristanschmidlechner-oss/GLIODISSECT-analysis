function s=pone_diag_state()
s=struct('attempts',0,'completed',0,'criterion',false,'bound',false,'fallback',uint16(0),'invalid',false,'caught',false,'skipped',false,'lastBaseline',NaN,'lastStep',NaN,'ratioBaseline',NaN,'ratioStep',NaN);
end
