function [fastq, vquest, QualityScoreAlign] = postVquestQfilter(fastqFile, vquestFile, vquest4file, vquest8file, start1,stop1,q1,start2,stop2,q2,name)
    
    delimiter = '\t';
    startRow = 1;
   	
   	% numbers of columns in each type of file
	n2 = 18;
	n4 = 18;
	n8 = 130;

	% format strings
	fmt2 = repmat('%s', 1, n2);
	fmt4 = repmat('%s', 1, n4);
	fmt8 = repmat('%s', 1, n8);
   	
    %formatSpec_2 = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
	%formatSpec_4 = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
	%formatSpec_8 = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

    disp('Reading Vquest File')
    disp(vquestFile)
    
	% read in header first, then remainder of lines
	tic; fid = fopen(vquestFile); h2 = fgetl(fid); C2 = textscan(fid, fmt2, 'delimiter', '\t'); fclose(fid); toc
	tic; fid = fopen(vquest4file); h4 = fgetl(fid); C4 = textscan(fid, fmt4, 'delimiter', '\t'); fclose(fid); toc
	tic; fid = fopen(vquest8file); h8 = fgetl(fid); C8 = textscan(fid, fmt8, 'delimiter', '\t'); fclose(fid); toc

    %fileID = fopen(vquestFile,'r');
    %fileID4 = fopen(vquest4file,'r');
    %fileID8 = fopen(vquest8file,'r');

    %dataArray = textscan(fileID, fmt2, 'Delimiter', delimiter, 'ReturnOnError', false);
	%dataArray4 = textscan(fileID4, fmt4, 'Delimiter', delimiter, 'ReturnOnError', false);
	%dataArray8 = textscan(fileID8, fmt8, 'Delimiter', delimiter, 'ReturnOnError', false);
	
    %fclose(fileID);

	% convert to cell matrix
	vquest = [C2{1:end}];
	vquest4 = [C4{1:end}];
	vquest8 = [C8{1:end}];

    %vquest = [dataArray{1:end}];
    %vquest4 = [dataArray4{1:end}];
    %vquest8 = [dataArray8{1:end}];
    %header2 = vquest(1,:);
    %header4 = vquest4(1,:);
    %header8 = vquest8(1,:);
    
    %vquest = vquest(2:end,:);
    %vquest4 = vquest4(2:end,:);
    %vquest8 = vquest8(2:end,:);
    
    % split header line into cell array of strings
	h2 = regexp(h2, '\t', 'split');
	h4 = regexp(h4, '\t', 'split');
	h8 = regexp(h8, '\t', 'split');

	% for some reason, headers were adding a field; remove
	header2 = h2(1:(end - 1));
	header4 = h4(1:(end - 1));
	header8 = h8(1:(end - 1));
    
    disp('Number of Vquest Sequences')
    disp(length(vquest))

    clearvars filename delimiter formatSpec fileID dataArray ans dataArray4 dataArray8;
    a = {};
    b = {};
    
    for k = 1:length(vquest)
		tmp = regexp(vquest{k,2},'^(\w+\_)','tokens');
		a(k,1) = tmp{1};
		
        b{k,1} =  regexprep(vquest{k,7},'\.','');
        vquest{k,2} = regexprep(vquest{k,2},'^\w+\_','');
        vquest{k,2} = regexprep(vquest{k,2},'000000000-','');
        vquest{k,2} = regexprep(vquest{k,2},'\_\|$','');
    end
    
    vquest = [vquest a b];
    
    disp('Reading FASTQ File')
    disp(fastqFile);
    fastq = fastqread(fastqFile);
    disp('Number of FASTQ Sequences')
    disp(length(fastq))
    
    for i = 1:length(fastq)
        fastq(i).Header = regexprep(fastq(i).Header,'000000000-','');
        fastq(i).RevComp = seqrcomplement(fastq(i).Sequence);
        fastq(i).RevQual = flip(fastq(i).Quality);
    end
	
	productiveIndicies = find(strcmp(vquest(:,3),'productive'));
    vquestProductive  = vquest(productiveIndicies,:);
    vquestProductive4 = vquest4(productiveIndicies,:);
    vquestProductive8 = vquest8(productiveIndicies,:);
    
    
    nv = length(vquestProductive);
    
    disp('Number of productive sequences')
    disp(nv)
    
    disp('Mapping Alignments')
    
    alignments = cell(nv,1);
    quality = cell(nv,1);
    qualityAlign = cell(nv,1);
    fwdscore = nan(nv,1);
    revscore = nan(nv,1);
    fwdalign = cell(nv,1);
    revalign = cell(nv,1);
    
    [~,ia,ib] = intersect({fastq.Header},vquestProductive(:,2));
    
    for j = 1:nv %change to nv
        [fwdscore(j),fwdalign{j}] = nwalign(vquestProductive{ib(j),end},fastq(ia(j)).Sequence,'GAPOPEN', 20, 'EXTENDGAP', 2);
        [revscore(j),revalign{j}] = nwalign(vquestProductive{ib(j),end},fastq(ia(j)).RevComp,'GAPOPEN', 20, 'EXTENDGAP', 2);
    end
    
    disp('Adjusting quality scores to alignment')
    [score,dmax] = max([fwdscore revscore],[],2);
    alignments(dmax == 1) = fwdalign(dmax == 1);
    alignments(dmax == 2) = revalign(dmax == 2);
    
    y = alignments; % change to alignments when above is changed to nv
    
    for j = 1:nv % change to nv
        tmp = fastq(ia(j)).Quality;
        a = alignments{j};
        tmp(strfind(a(1,:),'-')) = ' ';
        
        quality(j) = cellstr(tmp);
    end
    
    
    QualityNoGap = regexprep(quality,' ','');
    
    
    for j = 1:nv %change to nv
        a = vquestProductive{ib(j),7};
        b = QualityNoGap{j};
        
        n = length(a);
        q = setdiff(1:n,strfind(a,'.'));
        c = repmat(' ',1,n);
        c(q) = b;
        qualityAlign(j) = cellstr(c);
    end
  
    qualityScoreAlign = cellfun(@(x) double(x)-33, qualityAlign, 'UniformOutput',false);
    maxLen = max(cell2mat(cellfun(@(x) length(x), qualityScoreAlign, 'UniformOutput',false)));
    
    x = cellfun(@(x) [x (-1*ones(1,maxLen-length(x)))],qualityScoreAlign, 'UniformOutput', false);
    
    QualityScoreAlign = cell2mat(x);
    QualityScoreAlign(QualityScoreAlign == -1) = NaN;
    baseviz_r2(QualityScoreAlign);
    print(strcat(name,'qualityScoreDistribution_raw.eps'),'-depsc', '-noui');
    baseqc_r2(QualityScoreAlign);
    print(strcat(name,'qualityHistogram_raw.eps'),'-depsc', '-noui');
    
    ikeep = qFilter( QualityScoreAlign, name, start1,stop1,q1,start2,size(QualityScoreAlign,2),q2);
    n = length(ikeep);
    
    
    % add header to filtered data cell matrices
    vf2_iso = vquestProductive(ikeep,19);
    seqid = strcat(vf2_iso,vquestProductive(ikeep,2));
    
	vf2 = [header2; vquestProductive(ikeep, 1) seqid vquestProductive(ikeep, 3:18)];
	vf4 = [header4; vquestProductive4(ikeep, :)];
	vf8 = [header8; vquestProductive8(ikeep, :)];
	

	% output file names
	o2 = regexp(vquestFile,'.*\/(.*?\/.*?$)','tokens');
	o4 = regexp(vquest4file,'.*\/(.*?\/.*?$)','tokens');
	o8 = regexp(vquest8file,'.*\/(.*?\/.*?$)','tokens');
	
    o2 = strcat(pwd,'/','postVquestFiltered/',o2{1});
    o4 = strcat(pwd,'/','postVquestFiltered/',o4{1});
    o8 = strcat(pwd,'/','postVquestFiltered/',o8{1});

	folder = regexp(o2, '(.*)\/','tokens');
    folder = folder{1,1}{1};
	mkdir(folder{1});

	% format strings for output
	fmto2 = [repmat('%s\t', 1, n2-1) '%s\n'];
	fmto4 = [repmat('%s\t', 1, n4-1) '%s\n'];
	fmto8 = [repmat('%s\t', 1, n8-1) '%s\n'];

	% open file handles
	fid2 = fopen(char(o2), 'w');
	fid4 = fopen(char(o4), 'w');
	fid8 = fopen(char(o8), 'w');

	% write output
	for i = 1:(n + 1)
		fprintf(fid2, fmto2, vf2{i, :});
		fprintf(fid4, fmto4, vf4{i, :});
		fprintf(fid8, fmto8, vf8{i, :});
	end
	
	% close output files
	fclose(fid2);
	fclose(fid4);
	fclose(fid8);
    
    
    
    
    
    
    
    
    
    
%    outputFile = regexp(vquestFile,'.*\/(.*?\/.*?$)','tokens');
%    outputFile4 = regexp(vquest4file,'.*\/(.*?\/.*?$)','tokens');
%    outputFile8 = regexp(vquest8file,'.*\/(.*?\/.*?$)','tokens');
    
%    outputFile = strcat(pwd,'/','postVquestFiltered/',outputFile{1});
%    outputFile4 = strcat(pwd,'/','postVquestFiltered/',outputFile4{1});
%    outputFile8 = strcat(pwd,'/','postVquestFiltered/',outputFile8{1});
    

    
%	outputFileID = fopen(char(outputFile),'w');
%	outputFileID4 = fopen(char(outputFile4),'w');
%	outputFileID8 = fopen(char(outputFile8),'w');

%    Vq = vquestProductive(ikeep,:);
%    Vq4 = vquestProductive4(ikeep,:);
%    Vq8 = vquestProductive8(ikeep,:);
        
%    fprintf(outputFileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',header2{:});
%    fprintf(outputFileID4,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',header4{:});
%    fprintf(outputFileID8,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',header8{:});
    
%    [nrows,ncols] = size(Vq);
%    for row = 1:nrows
%    	fprintf(outputFileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',Vq{row,1:end-1});
%    	fprintf(outputFileID4,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',Vq4{row,:});
%    	fprintf(outputFileID8,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',Vq8{row,:});
%    end
    
    exit;
end



