function [ chans ] = ReadBioSemiTriggerChannel( trigChan )
% This function recieves the BioSemi status channel (the last channel in the data matrix read from
% the BDF file), and extracts the individual control channels multiplexed within it. See the biosemi
% website (http://www.biosemi.com/faq/trigger_signals.htm) for documentation of this channel
% structure. 
% 
% This function takes into account the possibility that the status message will change in the middle
% the sample, resulting in a spurrious intermediate sample where only some of the bits have been 
% updated. 
% A trigger is considered valid if it is non-zero, it is different from the previous value, and 
% one of the following conditions is fulfilled:
% 1. The next value is idential (example: third value in [0 0 10 10 10 0 0])
% 2. The previous and next values are equal (example: third value in [0 0 10 0]).
% 3. The previous value is also different from the one before it, but the one before that and 
% the next are equal (example: third value in [0 8 10 0]). 
%
% Written by Edden Gerber, June 2013
% 
%
% CHANGE HISTORY: 
%
% Oct.7 2014, Edden: Added the third condition listed in the description
% (ok3 in the code). 
% 

% channel bit ranges
bTRIG = 1:8;
bTRIG_EXTENDED = 1:16;
bSTART_EPOCH = 17;
bSPEED = [18:20 22];
bCMS_WITHIN_RANGE = 21;
bLOW_BATT = 23;
bACTIVE_TWO_MK2 = 24;

% Triggers
v = ExtractBits(trigChan,bTRIG);
val_changed = [v(1)~=0 ; diff(v)~=0];
prev_val_changed = [false ; val_changed(1:end-1)];
val_nonzero = v~=0;
next_val_is_same = [~val_changed(2:end) ; true];
neighbors_are_equal = [ true; v(3:end)-v(1:end-2) == 0 ; true];
next_and_one_before_prev_are_equal = [true ; true ; v(4:end)-v(1:end-3) == 0 ; true];
ok1 = val_nonzero & val_changed & next_val_is_same; % change to a stable non-zero value
ok2 = val_nonzero & val_changed & neighbors_are_equal; % single-sample change and return to previous value
ok3 = val_nonzero & val_changed & prev_val_changed & next_and_one_before_prev_are_equal; % a single-sample change which is spread across two samples, then returning to the original value
ok = ok1 | ok2 | ok3;
v(~ok) = 0;
chans.Triggers = v;
% Previous version: 
% change_loc = find(diff([0 ; v]) > 0 & diff([v ; 0]) <= 0); % first condition = value increased; second condition = value does not increase in the next sample.
% t = zeros(size(v));
% t(change_loc) = v(change_loc+1);
% chans.Triggers = t;

% Extended triggers
v = ExtractBits(trigChan,bTRIG_EXTENDED);
val_changed = [v(1)~=0 ; diff(v)~=0];
val_nonzero = v~=0;
next_val_is_same = [~val_changed(2:end) ; true];
neighbors_are_equal = [ true; v(3:end)-v(1:end-2) == 0 ; true];
ok1 = val_nonzero & val_changed & next_val_is_same; % change to a stable non-zero value
ok2 = val_nonzero & val_changed & neighbors_are_equal; % single-sample change and return to previous value
ok3 = val_nonzero & val_changed & prev_val_changed & next_and_one_before_prev_are_equal; % a single-sample change which is spread across two samples, then returning to the original value
ok = ok1 | ok2 | ok3;
v(~ok) = 0;
chans.ExtTriggers = v;
% Previous version: 
% change_loc = find(diff([0 ; v]) > 0 & diff([v ; 0]) <= 0); % first condition = value increased; second condition = value does not increase in the next sample.
% t = zeros(size(v));
% t(change_loc) = v(change_loc+1);
% chans.ExtTriggers = t;

% Start epoch
v = ExtractBits(trigChan,bSTART_EPOCH);
change_loc = find(diff([0 ; v]) ~= 0 & diff([v ; 0]) == 0);
t = zeros(size(v));
t(change_loc) = v(change_loc+1);
chans.StartEpoch = t;

% Speed
v = ExtractBits(trigChan,bSPEED);
change_loc = find(diff([0 ; v]) ~= 0 & diff([v ; 0]) == 0);
t = zeros(size(v));
t(change_loc) = v(change_loc+1);
chans.Speed = t;

% CMS is within Range
v = ExtractBits(trigChan,bCMS_WITHIN_RANGE);
change_loc = find(diff([0 ; v]) ~= 0 & diff([v ; 0]) == 0);
t = zeros(size(v));
t(change_loc) = v(change_loc+1);
chans.CmsInRange = t;

% Low battery
v = ExtractBits(trigChan,bLOW_BATT);
change_loc = find(diff([0 ; v]) ~= 0 & diff([v ; 0]) == 0);
t = zeros(size(v));
t(change_loc) = v(change_loc+1);
chans.LowBattery = t;

% ActiveTwo MK2
v = ExtractBits(trigChan,bACTIVE_TWO_MK2);
change_loc = find(diff([0 ; v]) ~= 0 & diff([v ; 0]) == 0);
t = zeros(size(v));
t(change_loc) = v(change_loc+1);
chans.ActiveTwoMk2 = t;

end



function [ out ] = ExtractBits( in, bits)
% This function extracts the specified bits from each number in a vector. 
% 
% in - the input vector
% bits - the bits to be extracted. The least significant bit (LSB) is 1. 
% For example, the decimal 10 is '1010' in binary form, and so 
% ExtractBits(10,[1 2 4]) will extract and concatenate the first, second and fourth bits (counting 
% from the end), resulting in '110' which is 6 in decimal (this is the result of subtracting the 3rd
% bit, correcponding to 4, from 10). 
%
% Written by Edden Gerber, June 2013

bits = sort(bits,'ascend');
out = zeros(size(in));
for i=1:length(bits)
    b = floor(mod(in,2^(bits(i)))/2^(bits(i)-1));
    b = b .* 2^(i-1);
    out = out + b;
end
end
