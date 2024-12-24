function visual_source(s_real, Cortex, thresh)
if nargin<3
    thresh = 0.1;
end
figure('color','k')
try
hp=patch(gca,'vertices',Cortex.vert,'faces',Cortex.face,...
    'FaceColor',[0.75 0.75 0.75],'edgecolor','none','facelighting','gouraud'...
    ,'specularstrength',0.2,'ambientstrength',0.5,'diffusestrength',0.5,...
    'BackfaceLighting','lit','AmbientStrength',0.5,'SpecularExponent',1,...
    'SpecularColorReflectance',0.5,'EdgeLighting','gouraud','facealpha',0.5);
catch e
hp=patch(gca,'vertices',Cortex.Vertices,'faces',Cortex.Faces,...
    'FaceColor',[0.75 0.75 0.75],'edgecolor','none','facelighting','gouraud'...
    ,'specularstrength',0.2,'ambientstrength',0.5,'diffusestrength',0.5,...
    'BackfaceLighting','lit','AmbientStrength',0.5,'SpecularExponent',1,...
    'SpecularColorReflectance',0.5,'EdgeLighting','gouraud','facealpha',0.5);    
end
material dull
% camlight('headlight','infinite');
set(gca,'color','k','Xcolor',[0 0 0],'ycolor',[0 0 0],'zcolor',[0 0 0],'CameraPosition', [0 0 0],'CameraViewAngle',6);
view([ -90 90 ])
axis equal
try
light('Position',[-100,0,-100]*mean(sum(Cortex.vert.^2,2)));light('Position',[100,0,100]*mean(sum(Cortex.vert.^2,2)));
catch e
light('Position',[-100,0,-100]*mean(sum(Cortex.Vertices.^2,2)));light('Position',[100,0,100]*mean(sum(Cortex.Vertices.^2,2)));    
end
colorMap = jet;
cmax = max(s_real);cmin = min(0,min(s_real)); 

tlen = size(colorMap,1);
if min(s_real)>=0
    colorMap(1:floor(tlen*thresh),:) = repmat([.75 .75 .75],floor(tlen*thresh),1);
else
    cmax = max(abs(s_real)); cmin = -cmax;
    colorMap(floor(tlen*(0.5-thresh/2))+1:floor(tlen*(thresh/2+0.5)),:) = repmat([.75 .75 .75],(floor(tlen*(thresh/2+0.5))-floor(tlen*(0.5-thresh/2))),1);
end
cdata=colorMap(floor(min((s_real-cmin)/(cmax-cmin),1)*(length(colorMap)-1))+1,:);
colormap(colorMap);
set(hp,'AlphaDataMapping','none', 'FaceVertexCData',cdata, 'facecolor', 'interp');