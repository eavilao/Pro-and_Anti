function [calData] = simpleEyeCalibration(calibValues,eyeX,eyeY,eyeFs)
%calibrates the eye data in x and y dimensions separately. Uses the offset
%and mag corrections given in calibValues. Returns a calData structure with
%fields: eyeXcal,eyeYcal,eyeFscal,eyeTcal



eyeX = (eyeX- calibValues.xOffset)./calibValues.xMagRatio;
eyeY = (eyeY- calibValues.yOffset)./calibValues.yMagRatio;

% now downsample to 1K
eyeXds = resample(eyeX,1000,round(eyeFs));
eyeYds = resample(eyeY,1000,round(eyeFs));

calData.eyeXcal = eyeXds;
calData.eyeYcal = eyeYds;
calData.eyeFsCal = 1000;
calData.eyeTcal = linspace(0,length(eyeXds)/1000,length(eyeXds));
end