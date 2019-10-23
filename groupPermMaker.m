function permSign = groupPermMaker(numSubj)
%numSubj must be at least 10
rng(1)
signIdx = [-1 1];
parfor p = 1:10000
    for s = 1:numSubj
        permSign(s,p) = signIdx(randperm(2,1));
    end
end
clear p s signIdx
permSign = unique(permSign','rows')';
permSign = permSign(:,randperm(size(permSign,2)));
[~,loc] = ismember(ones(1,size(permSign,1)),permSign','rows');
if loc ~=0
    permSign(:,loc) = [];
end
permSign = permSign(:,1:1023);
