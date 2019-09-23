function save_figures()
cd('plots')
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
n_fig = length(FigList);

for k = 2 : n_fig
      FigHandle = figure(k);
      FigName   = sprintf('figure_%d', k);
      savefig(FigHandle, [FigName, '.fig']);
      print(FigHandle,[FigName, '.eps'],'-depsc')
      
end

cd('..')
end