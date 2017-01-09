

query = 'select lineageID, population, count(*) from sequence group by lineageID, population order by lineageID, count(*)';
conn = igseq_db(db);
a = get(fetch(exec(conn,query)),'data');
close(conn);

xd = find(strcmp(a(:,2),'M3S9'));
xd = a(xd,:);

yd = find(strcmp(a(:,2),'M3S2'));
yd = a(yd,:);

[lineages,ia,ib] = intersect(cell2mat(yd(:,1)), cell2mat(xd(:,1)));

%sizes = cell2mat([M3S2_06v6(:,3)' M3S9_06v4(:,3)']);

x = zeros(length(lineages),1);
y = zeros(length(lineages),1);

for i = 1:length(lineages)
    x(i) = yd{ia(i),3};
    y(i) = xd{ib(i),3};
end




db = '/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/all_test_2016_09_20_16_08_05.db';

query = 'select lineageID, population, count(*) from sequence group by lineageID, population order by lineageID, count(*)';
conn = igseq_db(db);
a = get(fetch(exec(conn,query)),'data');
close(conn);



xd = find(strcmp(a(:,2),'M3S9'));
xd = a(xd,:);

yd = find(strcmp(a(:,2),'M3S2'));
yd = a(yd,:);

[lineages,ia,ib] = intersect(cell2mat(yd(:,1)), cell2mat(xd(:,1)));

%sizes = cell2mat([M3S2_06v6(:,3)' M3S9_06v4(:,3)']);

x = zeros(length(lineages),1);
y = zeros(length(lineages),1);

for i = 1:length(lineages)
    x(i) = yd{ia(i),3};
    y(i) = xd{ib(i),3};
end

pops = igseq_upop(db);
overlap = zeros(length(pops));
for i = 1:length(pops)
    for j = (i+1):length(pops)
        xd = find(strcmp(a(:,2),pops{i}));
        xd = a(xd,:);
        
        yd = find(strcmp(a(:,2),pops{j}));
        yd = a(yd,:);
        
        [lineages,ia,ib] = intersect(cell2mat(xd(:,1)), cell2mat(yd(:,1)));

        overlap(i,j) = length(lineages);
    end
end




xd = find(strcmp(a(:,2),'M3S7'));
xd = a(xd,:);

yd = find(strcmp(a(:,2),'M3S14'));
yd = a(yd,:);

[lineages,ia,ib] = intersect(cell2mat(yd(:,1)), cell2mat(xd(:,1)));

%sizes = cell2mat([M3S2_06v6(:,3)' M3S9_06v4(:,3)']);

x = zeros(length(lineages),1);
y = zeros(length(lineages),1);

for i = 1:length(lineages)
    x(i) = yd{ia(i),3};
    y(i) = xd{ib(i),3};
end




db = '/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/41_HepB_all_2016_09_28_14_23_14.db';

query = 'select lineageID, population, count(*) from sequence group by lineageID, population order by lineageID, count(*)';
conn = igseq_db(db);
a = get(fetch(exec(conn,query)),'data');
close(conn);



xd = find(strcmp(a(:,2),'M3S7'));
xd = a(xd,:);

yd = find(strcmp(a(:,2),'M3S14'));
yd = a(yd,:);

[lineages,ia,ib] = intersect(cell2mat(yd(:,1)), cell2mat(xd(:,1)));

%sizes = cell2mat([M3S2_06v6(:,3)' M3S9_06v4(:,3)']);

x = zeros(length(lineages),1);
y = zeros(length(lineages),1);

for i = 1:length(lineages)
    x(i) = yd{ia(i),3};
    y(i) = xd{ib(i),3};
end