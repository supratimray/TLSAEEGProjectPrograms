% code by Supratim; Modified by MD
function makeBox(h,xRange,yRange,colorName,lineWidth,lineStyle,boxType)

if ~exist('colorName','var') || isempty(colorName);           colorName=[0 0 0];              end
if ~exist('lineWidth','var') || isempty(lineWidth);           lineWidth=0.5;                  end
if ~exist('lineStyle','var') || isempty(lineStyle);           lineStyle='-';                  end

% H: Horizontal lines; V: Vertical lines; B: Box; M: Upper and Right margins (to complete a box for an axis)
if ~exist('boxType','var') || isempty(boxType);           boxType='H';                  end 

N=100;
xLine = xRange(1):diff(xRange)/N:xRange(2);
yLine = yRange(1):diff(yRange)/N:yRange(2);

hold(h,'on');
    switch boxType
        case 'H'
            plot(h,xLine,yRange(1)+zeros(1,N+1),'color',colorName,'lineWidth',lineWidth,'lineStyle',lineStyle);
            plot(h,xLine,yRange(2)+zeros(1,N+1),'color',colorName,'lineWidth',lineWidth,'lineStyle',lineStyle);
        case 'V'
            plot(h,xRange(1)+zeros(1,N+1),yLine,'color',colorName,'lineWidth',lineWidth,'lineStyle',lineStyle);
            plot(h,xRange(2)+zeros(1,N+1),yLine,'color',colorName,'lineWidth',lineWidth,'lineStyle',lineStyle);
        case 'B'
            plot(h,xLine,yRange(1)+zeros(1,N+1),'color',colorName,'lineWidth',lineWidth,'lineStyle',lineStyle);
            plot(h,xLine,yRange(2)+zeros(1,N+1),'color',colorName,'lineWidth',lineWidth,'lineStyle',lineStyle);
            plot(h,xRange(1)+zeros(1,N+1),yLine,'color',colorName,'lineWidth',lineWidth,'lineStyle',lineStyle);
            plot(h,xRange(2)+zeros(1,N+1),yLine,'color',colorName,'lineWidth',lineWidth,'lineStyle',lineStyle);  
        case 'M'
            plot(h,xLine,yRange(2)+zeros(1,N+1),'color',colorName,'lineWidth',lineWidth,'lineStyle',lineStyle);
            plot(h,xRange(2)+zeros(1,N+1),yLine,'color',colorName,'lineWidth',lineWidth,'lineStyle',lineStyle);
    end
end