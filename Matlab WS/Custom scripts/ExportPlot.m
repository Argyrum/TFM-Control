function ExportPlot(plotName, scaleX, scaleY, figureHandle)
%EXPORTPLOT (plotName, scaleX=1, scaleY=1, figureHandle=gcf)
%   Exports figureHandle to ./Plot out/ as .png
arguments (Input)
    plotName char
    scaleX (1,1) {mustBeNumeric} = 1
    scaleY (1,1) {mustBeNumeric} = 1
    figureHandle = gcf
end

theme(figureHandle,'light');
exportgraphics(figureHandle, ['./Plot out/', plotName , '.png'], Units="centimeters", ...
     Width=13.5*scaleX, Height=9.5*scaleY, Resolution=300)
end