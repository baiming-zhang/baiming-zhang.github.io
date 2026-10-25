'use strict';
(() => {
  const baseDraw=draw;let overlay;
  overlay=new FEBoundaryOverlay({canvas,gl,redraw:()=>draw()});
  draw=function(){
    if(gl&&program)gl.useProgram(program);baseDraw();if(!program||buffers.length!==models.length)return;
    const width=Math.floor(canvas.width/2),height=Math.floor(canvas.height/2),ca=Math.cos(az),sa=Math.sin(az),ce=Math.cos(el),se=Math.sin(el);
    const frames=[[0,'rodv2_0006-coarse','cal_peak','粗 Cal'],[1,'rodv2_0006','cal_peak','细 Cal'],[2,'rodv2_0006','cal_peak','细 Cal（显示平均）'],[3,'rodv2_0006','fen_peak','GPU FEn']].map(([panel,id,sample,title])=>({id,sample,title,center:common.center,right:[ca,sa,0],up:[-sa*se,ca*se,ce],depth:[sa*ce,-ca*ce,se],pan,scale:common.radius*.70/zoom,radius:common.radius,deform:1,viewport:{x:(panel%2)*width,y:panel<2?height:0,w:width,h:height}}));
    overlay.render(frames,program);
  };draw();
})();
