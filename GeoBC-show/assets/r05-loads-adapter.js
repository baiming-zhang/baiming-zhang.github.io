/* Load overlay for the retained twenty-model inspector, at its exact camera/scale. */
(()=>{'use strict';
const applied=new RodAppliedLoads(gl),data=new Map(),pending=new Set(),baseDraw=draw;let enabled=true;
const button=document.createElement('button');button.type='button';button.textContent='载荷';button.className='rod-load-button';button.setAttribute('aria-pressed','true');document.querySelector('.toolbar').append(button);
const legend=document.createElement('span');legend.className='rod-load-legend';legend.textContent='橙：销孔压力 · 紫红：惯性体载荷 · 箭长仅示意';document.querySelector('.foot').before(legend);
button.onclick=()=>{enabled=!enabled;button.setAttribute('aria-pressed',String(enabled));legend.hidden=!enabled;draw()};
draw=function(){gl.useProgram(program);baseDraw();if(!enabled||!active.length)return;const layout=grid(active.length),cw=Math.floor(canvas.width/layout.cols),ch=Math.floor(canvas.height/layout.rows),ca=Math.cos(az),sa=Math.sin(az),ce=Math.cos(el),se=Math.sin(el),radius=Math.max(...active.map(x=>x.geom.radius));
  active.forEach((entry,i)=>{let d=data.get(entry.caseId);if(!d){if(!pending.has(entry.caseId)){pending.add(entry.caseId);applied.load(entry.caseId).then(x=>{data.set(entry.caseId,x);draw()}).catch(e=>note('载荷标记：'+e.message))}return;}
    applied.draw({center:entry.geom.center,right:[ca,sa,0],up:[-sa*se,ca*se,ce],depth:[sa*ce,-ca*ce,se],pan,scale:radius*Math.max(1.2,1.2/(cw/ch))/zoom,radius,depthScale:radius*4,viewport:{x:(i%layout.cols)*cw,y:canvas.height-(Math.floor(i/layout.cols)+1)*ch,w:cw,h:ch}},d,'FE')
  })
};draw();
})();
