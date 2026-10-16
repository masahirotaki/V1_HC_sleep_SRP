function S = day5_test_sessions()
% DAY5_TEST_SESSIONS  Day-5 test days (ABCD, A_CD, E_CD back to back) per protocol.
%
% protocol 150 -> GAVNIK_<cond>, 250 -> GAVNIK250_<cond>. week = training week in
% which the animal ran that protocol (paired mice ran 150 and 250 in different weeks,
% see manuscript Fig. 5 legend). All recordings of a day share one sorting, so units
% are matched across conditions within a day (not across days).
%   E_CD: E is a novel orientation (>= 30 deg from any used before); grey = orientation 0.

% M00013 / M00014 are NPX1.0; M00013 ran TRAIN (150 ms + grey gaps) the week before.
% Not listed: M00071's later re-test days (20251202, 20251211) and GAVNIK200 (20251125).
C = {
% subject   date        protocol week
 'M00014', '20250321', 150, 1
 'M00013', '20250221', 150, 2
 'M00069', '20251003', 150, 1
 'M00071', '20251121', 150, 1
 'M00098', '20260512', 150, 1
 'M00088', '20260313', 150, 1
 'M00113', '20260716', 150, 1
 'M00087', '20260222', 150, 2
 'M00096', '20260703', 150, 2
 'M00100', '20260612', 150, 2
 'M00087', '20260217', 250, 1
 'M00096', '20260626', 250, 1
 'M00100', '20260606', 250, 1
 'M00088', '20260320', 250, 2
 'M00113', '20260724', 250, 2
};
S = struct('subject', C(:,1), 'date', C(:,2), 'protocol', C(:,3), 'week', C(:,4));
end
