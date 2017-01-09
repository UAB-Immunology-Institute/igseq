function baseviz_r2(seqs)

seqs(seqs == -1) = NaN;
seqs(seqs == 0) = NaN;

m = nanmean(seqs, 2);

[~, isort] = sort(m, 'descend');

[ns, n] = size(seqs);

pw = 1400;
ph = 800;
left = 50;
right = 100;
sp = 25;
idh = 20;
top = 20;
bot = 5;
fw = left + pw + right;
fh = bot + idh + sp + ph + top;

figure('position', [20 20 fw fh]);

pos1 = [left / fw, bot / fh, pw / fw, idh / fh];
pos2 = [left / fw, (bot + idh + sp) / fh, pw / fw, ph / fh];

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
seqs(isnan(seqs)) = 0;
imagesc(seqs(isort, :), [0 max(seqs(:))]);
colormap(jet1(256));
plot([1 1] * max(fr1),  ylim, 'w-', 'linewidth', 1.25);
plot([1 1] * max(cdr1), ylim, 'w-', 'linewidth', 1.25);
plot([1 1] * max(fr2),  ylim, 'w-', 'linewidth', 1.25);
plot([1 1] * max(cdr2), ylim, 'w-', 'linewidth', 1.25);
plot([1 1] * max(fr3),  ylim, 'w-', 'linewidth', 1.25);
set(gca,...
    'xlim', [.5 n + .5],...
    'ylim', [.5 ns + .5],...
    'ticklength', [0 0],...
    'fontname', 'arial',...
    'fontsize', 11);

slider = uicontrol(...
    'style', 'slider',...
    'visible', 'on',...
    'min', 1,...
    'max', ns - 400,...
    'value', 1,...
    'units', 'normalized',...
    'position', [(fw - 40) / fw, (bot + idh + sp) / fh, 20 / fw, ph / fh],...
    'sliderstep', [1 / ns, 400 / ns],...
    'callback', @cb_slider); %#ok<NASGU>

button = uicontrol(...
    'style', 'pushbutton',...
    'visible', 'on',...
    'units', 'normalized',...
    'position', [(fw - 80) / fw, bot / fh, 60 / fw, idh / fh],...
    'string', 'RESET',...
    'callback', @cb_reset); %#ok<NASGU>

    H.ax = gca;
    H.ns = ns;

    set(gcf, 'userdata', H);
    
return


function cb_slider(hcb, ~)
    h = get(gcf, 'userdata');
    k = get(hcb, 'value');
    set(h.ax, 'ylim', [k - .5 k + 400.5]);
return

function cb_reset(~, ~)
    h = get(gcf, 'userdata');
    set(h.ax, 'ylim', [.5 h.ns + .5]);
return




