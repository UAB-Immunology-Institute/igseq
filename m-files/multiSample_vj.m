fn = {...
    'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S14_Human_a_v6_2016_07_28_11_06_47.db';...
    'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S25_Human_a_v6_2016_07_28_13_25_10.db';...
    'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S26_Human_a_v6_2016_07_28_13_39_22.db';...
    'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S15_Human_a_v6_2016_07_28_12_23_16.db';...
    'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S2_Human_a_v6_2016_07_28_09_36_12.db';...
    'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M3S5_Human_a_v6_2016_07_28_10_14_44.db';...
    'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/M6S40_Human_a_v6_2016_07_28_15_29_47.db'};

lab = regexprep(fn, 'Users/arosenberg/desktop/IGSEQ_DEV/regeneron/', '');

for i = 1:7
    temp = lab{i};
    k = min(strfind(temp, '_'));
    lab(i) = cellstr(temp(1:(k - 1)));
end

S = struct;
L = struct;
SG = struct;
LG = struct;
uv = {};
for i = 1:7
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

vjseq = zeros(7, nv, nj);
vjlin = zeros(7, nv, nj);
for i = 1:7
    tempseq = zeros(nv, nj);
    templin = zeros(nv, nj);
    [~, ~, idx] = intersect(S(i).vj(3:end, 1), uv);
    tempseq(idx, :) = cell2mat(S(i).vj(3:end, 2:end));
    templin(idx, :) = cell2mat(L(i).vj(3:end, 2:end));
%     [~, ~, idx] = intersect(SG(i).vj(3:end, 1), uv);
%     tempseq(idx, :) = cell2mat(SG(i).vj(3:end, 2:end));
%     templin(idx, :) = cell2mat(LG(i).vj(3:end, 2:end));
    vjseq(i, :, :) = 100 * tempseq / sum(tempseq(:));
    vjlin(i, :, :) = 100 * templin / sum(templin(:));
end

pix = 15;
pw = pix * nj;
ph = pix * nv;
bot = 30;
right = 70;
top = 40;
fw = (pw + right) * 7;
fh = bot + ph + top;
pos  = (0:6) * (pw + right) / fw;
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
for f = 1:7   
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
for f = 1:7   
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
    
    
vjseqprof = nan(nv * nj, 7);
vjlinprof = nan(nv * nj, 7);

for i = 1:7
    temp = squeeze(vjseq(i, :, :))';
    vjseqprof(:, i) = temp(:);
    temp = squeeze(vjlin(i, :, :))';
    vjlinprof(:, i) = temp(:);
end

jj = repmat(uj', nv, 1);
vv = repmat(uv', nj, 1);
vv = vv(:);







% isotypes

iso = nan(7, 3);
for i = 1:7
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
hd = nan(length(x), 7);
he = nan(length(x), 7);

for i = 1:7
    hd(:, i) = igseq_hillplot(fn{i}, x, 0, lab(i), {''}, {'r'});
    he(:, i) = igseq_hillplot(fn{i}, x, 1, lab(i), {''}, {'b'});
end

figure;
hold on
plot(x, he(:, 1:3), '-', 'linewidth', 1.5, 'color', [1.0 0.6 0.0]);
plot(x, he(:, 4:7), '-', 'linewidth', 1.5, 'color', [0.5 0.5 1.0]);





    
