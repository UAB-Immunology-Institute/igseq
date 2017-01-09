function baseqc_r2(seqs)



m2 = seqs;
m2(m2 == 0) = NaN;
m2(m2 == -1) = NaN;

[ns, n] = size(seqs);

edges = -.5:1:40.5;

h = histc(m2, edges);

h = h(1:(end - 1), :);

s = nansum(h);

z = zeros(n, 1);
for i = 1:ns
    k = find(seqs(i, :) == 0, 1);
    if ~isempty(k), z(k) = z(k) + 1; else z(n) = z(n) + 1; end
end   

pw = 1400;
ph = 250;
left = 50;
right = 125;
sp = 25;
gh = 90;
idh = 20;
top = 8;
bot = 5;
fw = left + pw + right;
fh = bot + idh + sp + ph + sp + gh + sp + gh + top;
ncol = 256;
cmap = jet(ncol);
cmap(1, :) = [.7 .7 .7];
cmax = max(log10(1 + h(:)));

figure('position', [20 20 fw fh]);

pos1 = [left / fw, bot / fh, pw / fw, idh / fh];
pos2 = [left / fw, (bot + idh + sp) / fh, pw / fw, ph / fh];
pos3 = [left / fw, (bot + idh + sp + ph + sp) / fh, pw / fw, gh / fh];
pos4 = [left / fw, (bot + idh + sp + ph + sp + gh + sp) / fh, pw / fw, gh / fh];

% V gene segments
fr1  = [.5    78.5];
cdr1 = [78.8  114.5];
fr2  = [114.5 165.5];
cdr2 = [165.5 195.5];
fr3  = [195.5 312.5];
cdr3 = [312.5 n + .5];

% V gene ideogram
subplot('position', pos1, 'nextplot', 'add');
fill(fr1([1 2 2 1 1]), [0 0 1 1 0], [.9 .6 .9], 'edgecolor', 'none');
fill(fr2([1 2 2 1 1]), [0 0 1 1 0], [.9 .6 .9], 'edgecolor', 'none');
fill(fr3([1 2 2 1 1]), [0 0 1 1 0], [.9 .6 .9], 'edgecolor', 'none');
fill(cdr1([1 2 2 1 1]), [0 0 1 1 0], [1 .8 .4], 'edgecolor', 'none');
fill(cdr2([1 2 2 1 1]), [0 0 1 1 0], [1 .8 .4], 'edgecolor', 'none');
fill(cdr3([1 2 2 1 1]), [0 0 1 1 0], [.7 .7 .7], 'edgecolor', 'none');
text(mean(fr1),  .5, 'FR1',        'fontname', 'arial', 'fontsize', 16, 'horizontalalignment', 'center', 'fontweight', 'bold');
text(mean(fr2),  .5, 'FR2',        'fontname', 'arial', 'fontsize', 16, 'horizontalalignment', 'center', 'fontweight', 'bold');
text(mean(fr3),  .5, 'FR3',        'fontname', 'arial', 'fontsize', 16, 'horizontalalignment', 'center', 'fontweight', 'bold');
text(mean(cdr1), .5, 'CDR1',       'fontname', 'arial', 'fontsize', 16, 'horizontalalignment', 'center', 'fontweight', 'bold');
text(mean(cdr2), .5, 'CDR2',       'fontname', 'arial', 'fontsize', 16, 'horizontalalignment', 'center', 'fontweight', 'bold');
text(mean(cdr3), .5, 'CDR3 + FR4', 'fontname', 'arial', 'fontsize', 16, 'horizontalalignment', 'center', 'fontweight', 'bold');
set(gca,...
    'ticklength', [0 0],...
    'xticklabel', {},...
    'yticklabel', {},...
    'xcolor', 'w',...
    'ycolor', 'w',...
    'box', 'off',...
    'ylim', [0 1],...
    'xlim', [.5 n + .5]);

% heat map
subplot('position', pos2, 'nextplot', 'add');
for i = 1:size(h, 1)
    for j = 1:n
        v = round((log10(1 + h(i, j)) / cmax) * (ncol - 1)) + 1;
        x = [-.5 .5 .5 -.5 -.5] + j;
        y = [-.5 -.5 .5 .5 -.5] + i;
        fill(x, y - 1, cmap(v, :), 'edgecolor', cmap(v, :), 'linewidth', .5);
    end
end
plot([1 1] * max(fr1),  [-.5 40.5], 'k-', 'linewidth', 1.25);
plot([1 1] * max(cdr1), [-.5 40.5], 'k-', 'linewidth', 1.25);
plot([1 1] * max(fr2),  [-.5 40.5], 'k-', 'linewidth', 1.25);
plot([1 1] * max(cdr2), [-.5 40.5], 'k-', 'linewidth', 1.25);
plot([1 1] * max(fr3),  [-.5 40.5], 'k-', 'linewidth', 1.25);
plot([.5 n + .5], [1 1] * 29.5, 'k-', 'linewidth', 1.25);
set(gca,...
    'ticklength', [0.000 0.000],...
    'tickdir', 'out',...
    'xlim', [.5 n + .5],...
    'ytick', 0:5:40,...
    'ylim', [-.5 40.5],...
    'box', 'off',...
    'fontname', 'arial',...
    'fontsize', 11);
text(1.005 * n, mean(ylim), {'Quality'; 'Score'}, 'fontname', 'arial', 'fontsize', 18);

% number of sequences bar plot
subplot('position', pos3, 'nextplot', 'add');
bar(1:n, s, 1, 'facecolor', [1 1 1] * .3, 'edgecolor', 'none');
set(gca,...
    'ticklength', [0 0],...
    'xticklabel', {},...
    'xlim', [.5 n + .5],...
    'box', 'off',...
    'fontname', 'arial',...
    'fontsize', 11);
set(gca, 'yticklabel', cellstr(num2str(get(gca, 'ytick')')));
plot([1 1] * max(fr1),  ylim, 'k-', 'linewidth', 1.25);
plot([1 1] * max(cdr1), ylim, 'k-', 'linewidth', 1.25);
plot([1 1] * max(fr2),  ylim, 'k-', 'linewidth', 1.25);
plot([1 1] * max(cdr2), ylim, 'k-', 'linewidth', 1.25);
plot([1 1] * max(fr3),  ylim, 'k-', 'linewidth', 1.25);
text(1.005 * n, mean(ylim), {'# Sequences'; 'per nt Position'}, 'fontname', 'arial', 'fontsize', 18);

% distributions of ending points
subplot('position', pos4, 'nextplot', 'add');
bar(1:n, z, 1, 'facecolor', [1 1 1] * .3, 'edgecolor', 'none');
set(gca,...
    'ticklength', [0 0],...
    'xticklabel', {},...
    'xlim', [.5 n + .5],...
    'box', 'off',...
    'fontname', 'arial',...
    'fontsize', 11);
set(gca, 'yticklabel', cellstr(num2str(get(gca, 'ytick')')));
plot([1 1] * max(fr1),  ylim, 'k-', 'linewidth', 1.25);
plot([1 1] * max(cdr1), ylim, 'k-', 'linewidth', 1.25);
plot([1 1] * max(fr2),  ylim, 'k-', 'linewidth', 1.25);
plot([1 1] * max(cdr2), ylim, 'k-', 'linewidth', 1.25);
plot([1 1] * max(fr3),  ylim, 'k-', 'linewidth', 1.25);
text(1.005 * n, mean(ylim), {'Distribution'; 'of Termini'}, 'fontname', 'arial', 'fontsize', 18);


return