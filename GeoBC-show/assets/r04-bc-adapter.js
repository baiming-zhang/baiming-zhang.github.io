'use strict';
(() => {
  const baseDraw=draw,overlay=new FEBoundaryOverlay({canvas,gl,redraw:()=>draw()});
  draw=function(){
    if(program)gl.useProgram(program);baseDraw();if(!program||!active.length)return;
    const layout=grid(active.length),cw=Math.floor(canvas.width/layout.cols),ch=Math.floor(canvas.height/layout.rows),ca=Math.cos(az),sa=Math.sin(az),ce=Math.cos(el),se=Math.sin(el),radius=Math.max(...active.map(x=>x.geom.radius)),maxU=Math.max(1e-12,...active.filter(x=>x.field).map(x=>x.sample.max_U_mm||0)),deformation=Math.min(8,.18*radius/maxU);
    const frames=active.map((entry,i)=>({id:entry.caseId,sample:scenario,captionOffset:86,showDeform:true,title:entry.caseId+' · '+(manifest.scenarios.find(x=>x.id===scenario)?.label||scenario),center:entry.geom.center,right:[ca,sa,0],up:[-sa*se,ca*se,ce],depth:[sa*ce,-ca*ce,se],pan,radius,scale:radius*Math.max(.78,1.05/(cw/ch))/zoom,deform:mode==='U'&&entry.field?deformation:0,viewport:{x:(i%layout.cols)*cw,y:canvas.height-(Math.floor(i/layout.cols)+1)*ch,w:cw,h:ch}}));
    overlay.render(frames,program);
  };draw();
})();
