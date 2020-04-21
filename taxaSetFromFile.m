function obj = taxaSetFromFile(file,cladeNames,nLdmk)

% Generate taxaSet from .csv file inputs:
%
% file = .csv file with .tif input file names, First Appearance Dates
%        (FAD), Last Appearance Dates (LAD) and clade numbers 
%
% row of input data, to load image 'Acanthostega.tif', with a FAD in the
% 3rd time bin, LAD in the 5th time bin, and belonging to clade 1
%
% |    TAXON     |  FAD  |  LAD  | CLADE |
% | Acanthostega |   3   |   5   |   1   |
%
%
% cladeNames = .csv file with list of clade names for plotting
%
% nLdmk = number of landmarks to inpu into EFA. This is standardised for
%           all taxa. Use a sensitivity analysis to obtain the optimal
%           value, which usually lies around 500 - 600, depending on image
%           quality/complexity.

input = importdata(file);
cNames = importdata(cladeNames);

names = input.textdata;

[N,~] = size(names);

FAD = input.data(:,1);
LAD = input.data(:,2);

clade = input.data(:,3);

Pwr = zeros(N,39);

taxa = taxonData.empty(0,1);

LoadBar(0,N);

for i = 1:N
    disp(names{i});
    n = numel(names{i}) + 1;
    taxa(i) = taxonData(names(i),FAD(i),LAD(i),clade(i),nLdmk);
    if ~isempty(taxa(i).pwr)
    	Pwr(i,:) = taxa(i).pwr;
    end
    LoadBar(i,N,n);
end

i = 1;
while i <= length(taxa)
    if isempty(taxa(i).harm)
        S = strcat('WARNING: ',names{i},' could not be landmarked and was omitted. Please check image file.');
        disp(S);
        taxa(i) = [];
        names(i) = [];
        Pwr(i,:) = [];
        i = i - 1;
    end
    i = i + 1;
end
    

sPwr = prctile(Pwr,[5 95]);

search = true;

for i = 1:39
    if sPwr(1,i) > 0.99 && search
        search = false;
        nHarm = i+1;
    end
end

harm = zeros(N,4*nHarm - 3);

for i = 1:length(taxa)
    harm(i,:) = taxa(i).harmRow(nHarm);
end


[c,s,~,~,e,m] = pca(harm);


obj = taxaSet(taxa,nHarm,s,c,e,m,cNames);

end


