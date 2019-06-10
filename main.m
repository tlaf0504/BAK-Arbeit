clear all
close all
clc

global ui;
ui = 1;

addpath(fullfile(pwd, 'entities'))

if ui
    main_window()
else
end

%tic
%Parser.parse('/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/cap_high_resolution', 'capacitor.geo')
%Parser.parse('/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/simple_cap', 'simple_cap.geo')
%toc


        