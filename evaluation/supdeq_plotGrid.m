%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function supdeq_plotGrid( samplingGrid )
%
% This function plots the spatial sampling grid on a sphere and on a 2D map
%
% Output:
% Figure
%
% Input:
% samplingGrid  - Spatial sampling grid (Q x 2 matrix), where the first 
%                 column holds the azimuth and the second the
%                 elevation (both in degree).
%                 Azimuth in degree (0=front, 90=left, 180=back, 270=right)
%                 (0 points to positive x-axis, 90 to positive y-axis)
%                 Elevations in degree (0=North Pole, 90=front, 180=South Pole)
%                 (0 points to positive z-axis, 180 to negative z-axis)
% plotMode      - 'default' or 'gray', 'gray' is independent on sample size
%                  and delivers no map plot
%                 Default: 'default'
%
% Dependencies: -
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

%%
function supdeq_plotGrid ( samplingGrid, plotMode )

if nargin < 2
    plotMode = 'default';
end

samplingGridRad(:,1) = samplingGrid(:,1) * pi / 180;
samplingGridRad(:,2) = samplingGrid(:,2) * pi / 180;
[Xm,Ym,Zm]=sph2cart(samplingGridRad(:,1),pi/2-samplingGridRad(:,2),1.01);


if strcmp(plotMode,'default')

    figure %Sphere plot
    colormap Gray;
    if size(Xm,1)>1500
        plot3(Xm,Ym,Zm,'marker','.','markerfacecolor','g','color','g','linestyle','none')
    else
        plot3(Xm,Ym,Zm,'marker','o','markerfacecolor','g','color','g','linestyle','none')
    end

    axis off;
    hold on;
    grid off;
    sphere;
    axis equal;
    rotate3d on;
    light;
    alpha(.8);
    lighting phong;
    camzoom(1.4);
    view(0,0)
    hold off;

    figure %Map plot
    plot(samplingGrid(:,1),samplingGrid(:,2),'.')
    set(gca, 'YDir','reverse')
    xlim([-10 369])
    ylim([-10 190])
    xlabel('Azimuth [°]');
    ylabel('Elevation [°]');

elseif strcmp(plotMode,'gray')
    
    %figure %Sphere plot
    colormap Gray;
    
    plot3(Xm,Ym,Zm,'marker','o','markerfacecolor','k','color','k','linestyle','none','markersize', 1.5*8)

    axis off; 
    grid off; 
    hold on; 
    sphere; 
    shading interp; 
    axis equal; 
    rotate3d off; 
    %alpha(0.75); 
    lighting gouraud; 
    hold off; 

end


end

