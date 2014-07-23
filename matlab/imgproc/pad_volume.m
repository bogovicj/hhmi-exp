function [volpad, origRng] = pad_volume(vol, amt, padvalue)
% volpad = pad_volume(vol,amt, padvalue)

if(~exist('padvalue','var') || isempty(padvalue))
   padvalue = 0; 
end

if(isscalar(amt))
   padVec =  repmat(amt, [1 3]);
else
   padVec = amt;
end

sz = size(vol);
szPad = sz + 2.*padVec;

volpad = padvalue.*ones(szPad);
volpad( padVec(1)+1: padVec(1) + sz(1), ...
        padVec(2)+1: padVec(2) + sz(2), ... 
        padVec(3)+1: padVec(3) + sz(3)) = vol;
  

origRng = [	amt+1; amt+sz(1); ...
			amt+1; amt+sz(2); ...
			amt+1; amt+sz(3)]; 
