function dif = dif_RCV(rcv)
%DIF_RCV Mean paired difference of the robust CV (RCV) between two conditions.
%
%   dif = dif_RCV(rcv) returns the mean, across subjects, of the difference
%   between the RCV of the first and the second column of rcv. It is the
%   statistic used in the bootstrap comparison of the examples
%   (see Examples/example_statistics.m).
%
% Inputs:
%   rcv: [subjects x 2] matrix with the RCV of each subject in two
%        conditions (column 1: first condition, column 2: second condition).
%
% Outputs:
%   dif: mean difference rcv(:,1) - rcv(:,2) (scalar).

dif = rcv(:,1) - rcv(:,2);
dif = mean(dif);
end
