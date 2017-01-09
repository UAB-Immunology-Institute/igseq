% fn = {...
%     'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S14_Human_a_v6_2016_07_28_11_06_47.db';...
%     'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S25_Human_a_v6_2016_07_28_13_25_10.db';...
%     'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S26_Human_a_v6_2016_07_28_13_39_22.db';...
%     'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S15_Human_a_v6_2016_07_28_12_23_16.db';...
%     'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S2_Human_a_v6_2016_07_28_09_36_12.db';...
%     'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S5_Human_a_v6_2016_07_28_10_14_44.db';...
%     'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M6S40_Human_a_v6_2016_07_28_15_29_47.db'};

fn = {...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M1S1_Human_b_v4_2016_08_23_13_48_23.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S16_Human_b_v4_2016_08_25_09_53_30.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M1S2_Human_b_v4_2016_08_23_14_37_08.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S11_Human_b_v4_2016_08_25_10_18_40.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M1S3_Human_b_v4_2016_08_23_15_42_25.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S9_Human_b_v4_2016_08_25_10_38_04.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M1S4_Human_b_v4_2016_08_23_16_31_03.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S19_Human_b_v4_2016_08_25_10_56_28.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M1S5_Human_b_v4_2016_08_24_09_33_04.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S20_Human_b_v4_2016_08_25_11_25_37.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S22_Human_b_v4_2016_08_25_15_56_45.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S24_Human_b_v4_2016_08_25_15_25_31.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S6_Human_b_v4_2016_08_25_14_34_47.db';...

'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S23_Human_b_v4_2016_08_25_15_04_51.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M1S6_Human_b_v4_2016_08_24_10_37_11.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S13_Human_b_v4_2016_08_25_12_03_02.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M1S7_Human_b_v4_2016_08_24_11_28_08.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S7_Human_b_v4_2016_08_25_12_26_56.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M1S8_Human_b_v4_2016_08_24_12_34_00.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S1_Human_b_v4_2016_08_25_13_01_13.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M1S9_Human_b_v4_2016_08_24_15_59_42.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S17_Human_b_v4_2016_08_25_13_23_38.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M1S10_Human_b_v4_2016_08_25_08_54_03.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S12_Human_b_v4_2016_08_25_14_05_13.db';...
'/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/M3S4_Human_b_v4_2016_08_25_16_45_44.db'};

n = length(fn);

lab = regexprep(fn, '/Volumes/G-DRIVE/Kobie/Regeneron_analysis/Human/Human_b_V4_PB/db/', '');

for i = 1:n
    temp = lab{i};
    k = min(strfind(temp, '_'));
    lab(i) = cellstr(temp(1:(k - 1)));
end

S = struct;
L = struct;
SG = struct;
LG = struct;
uv = {};
for i = 1:n
    S(i).vj = igseq_vjcount(fn{i}, 'sequence', 1);
    L(i).vj = igseq_vjcount(fn{i}, 'lineage', 1);
    SG(i).vj = igseq_vjcount_IGG(fn{i}, 'sequence', 1);
    LG(i).vj = igseq_vjcount_IGG(fn{i}, 'lineage', 1);
    uv = [uv; S(i).vj(3:end, 1)];
end
uv = unique(uv);
uj = S(1).vj(2, 2:end);

uv0 = regexprep(uv, 'IGHV', 'V');
uj0 = regexprep(uj, 'IGHJ', 'J');

nv = length(uv);
nj = length(uj);

vjseq = zeros(n, nv, nj);
vjlin = zeros(n, nv, nj);
for i = 1:n
    tempseq = zeros(nv, nj);
    templin = zeros(nv, nj);
%     [~, ~, idx] = intersect(S(i).vj(3:end, 1), uv);
%     tempseq(idx, :) = cell2mat(S(i).vj(3:end, 2:end));
%     templin(idx, :) = cell2mat(L(i).vj(3:end, 2:end));
     [~, ~, idx] = intersect(SG(i).vj(3:end, 1), uv);
     tempseq(idx, :) = cell2mat(SG(i).vj(3:end, 2:end));
     templin(idx, :) = cell2mat(LG(i).vj(3:end, 2:end));
    vjseq(i, :, :) = 100 * tempseq / sum(tempseq(:));
    vjlin(i, :, :) = 100 * templin / sum(templin(:));
end

pix = 15;
pw = pix * nj;
ph = pix * nv;
bot = 30;
right = 70;
top = 40;
fw = (pw + right) * n;
fh = bot + ph + top;
pos  = (0:n-1) * (pw + right) / fw;
pos1 = pos + (pw / fw);
pref = struct(...
    'ticklength', [0 0],...
    'xticklabel', {{}},...
    'yticklabel', {{}},...
    'xcolor', 'w',...
    'ycolor', 'w',...
    'box', 'off',...
    'ydir', 'reverse');
ncol = 1024;
cmap = parula(ncol);

cmaxlin = 7;
cmaxseq = 20;

figure('position', [20 20 fw fh]);
for f = 1:n   
    subplot('position', [pos1(f), bot / fh, right / fw, ph / fh]);
        for k = 1:nv, text(.1, k, uv0{k}, 'fontname', 'arial', 'fontsize', 10); end
        set(gca, pref, 'xlim', [0 1], 'ylim', [.5 nv + .5]);    
    subplot('position', [pos(f), 0, pw / fw, bot / fh]);
        for k = 1:nj, text(k, .15, uj0{k}, 'fontname', 'arial', 'fontsize', 10, 'rotation', 270); end
        set(gca, pref, 'xlim', [.5 nj + .5], 'ylim', [0 1]);  
    subplot('position', [pos(f), (bot + ph) / fh, pw / fw, top / fh]);
        text(.5, .5, lab{f}, 'fontname', 'arial', 'fontsize', 18, 'fontweight', 'bold', 'horizontalalignment', 'center');
        set(gca, pref, 'xlim', [0 1], 'ylim', [0 1]);
    subplot('position', [pos(f), bot / fh, pw / fw, ph / fh], 'nextplot', 'add');
        data = squeeze(vjseq(f, :, :));
        for i = 1:nv
            for j = 1:nj
                if data(i, j) > 0
                    v = round((data(i, j) / cmaxseq) * (ncol - 1)) + 1;
                    if v > ncol, v = ncol; end
                    c = cmap(v, :);
                else
                    c = [0 0 0];
                end
                x = [-.5 .5 .5 -.5 -.5] + j;
                y = [-.5 -.5 .5 .5 -.5] + i;
                fill(x, y, c, 'edgecolor', c);
            end
        end
        set(gca, pref, 'xlim', [.5 nj + .5], 'ylim', [.5 nv + .5]);
end
    
figure('position', [20 20 fw fh]);
for f = 1:n 
    subplot('position', [pos1(f), bot / fh, right / fw, ph / fh]);
        for k = 1:nv, text(.1, k, uv0{k}, 'fontname', 'arial', 'fontsize', 10); end
        set(gca, pref, 'xlim', [0 1], 'ylim', [.5 nv + .5]);    
    subplot('position', [pos(f), 0, pw / fw, bot / fh]);
        for k = 1:nj, text(k, .15, uj0{k}, 'fontname', 'arial', 'fontsize', 10, 'rotation', 270); end
        set(gca, pref, 'xlim', [.5 nj + .5], 'ylim', [0 1]);  
    subplot('position', [pos(f), (bot + ph) / fh, pw / fw, top / fh]);
        text(.5, .5, lab{f}, 'fontname', 'arial', 'fontsize', 18, 'fontweight', 'bold', 'horizontalalignment', 'center');
        set(gca, pref, 'xlim', [0 1], 'ylim', [0 1]);
    subplot('position', [pos(f), bot / fh, pw / fw, ph / fh], 'nextplot', 'add');
        data = squeeze(vjlin(f, :, :));
        for i = 1:nv
            for j = 1:nj
                if data(i, j) > 0
                    v = round((data(i, j) / cmaxlin) * (ncol - 1)) + 1;
                    if v > ncol, v = ncol; end
                    c = cmap(v, :);
                else
                    c = [0 0 0];
                end
                x = [-.5 .5 .5 -.5 -.5] + j;
                y = [-.5 -.5 .5 .5 -.5] + i;
                fill(x, y, c, 'edgecolor', c);
            end
        end
        set(gca, pref, 'xlim', [.5 nj + .5], 'ylim', [.5 nv + .5]);
end
    
    
vjseqprof = nan(nv * nj, n);
vjlinprof = nan(nv * nj, n);

for i = 1:n
    temp = squeeze(vjseq(i, :, :))';
    vjseqprof(:, i) = temp(:);
    temp = squeeze(vjlin(i, :, :))';
    vjlinprof(:, i) = temp(:);
end

jj = repmat(uj', nv, 1);
vv = repmat(uv', nj, 1);
vv = vv(:);


for i = 1:n
    igseq_runAll(fn{i});
end




% isotypes

iso = nan(n, 3);
for i = 1:n
    z = igseq_isotypes(fn{i});
    ty = z(1, 2:end);
    nm = cell2mat(z(2, 2:end));
    im = find(strncmp(ty, 'M', 1));
    ig = find(strncmp(ty, 'G', 1));
    ia = find(strncmp(ty, 'A', 1));
    iso(i, :) = [sum(nm(im)) sum(nm(ia)) sum(nm(ig))];
end
iso = 100 * iso ./ repmat(sum(iso, 2), 1, 3);



% hill

x = 0:.01:3;
hd = nan(length(x), n);
he = nan(length(x), n);

for i = 1:n
    hd(:, i) = igseq_hillplot(fn{i}, x, 0, lab(i), {''}, {'r'});
    he(:, i) = igseq_hillplot(fn{i}, x, 1, lab(i), {''}, {'b'});
end

figure;
hold on
plot(x, he(:, 1:13), '-', 'linewidth', 1.5, 'color', [1.0 0.6 0.0]);
plot(x, he(:, 14:n), '-', 'linewidth', 1.5, 'color', [0.5 0.5 1.0]);





    
