'use strict';
const canvas=document.getElementById('canvas'),gl=canvas.getContext('webgl2',{antialias:true});
let program,overlay,dataset,bodyTri,bodyCount,az=.65,el=.35,roll=0,zoom=1,pan=[0,0],drag=null;
const caseId=new URLSearchParams(location.search).get('case')==='r02-reference-patch'?'r02-reference-patch':'c01-smoke';
document.getElementById('case').value=caseId;document.getElementById('back').href=caseId==='c01-smoke'?'../C01.html':'../R02.html';
function shader(type,s){const q=gl.createShader(type);gl.shaderSource(q,s);gl.compileShader(q);if(!gl.getShaderParameter(q,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(q));return q;}
function init(){
  if(!gl)throw Error('需要 WebGL2 浏览器');program=gl.createProgram();
  gl.attachShader(program,shader(gl.VERTEX_SHADER,`#version 300 es
in vec3 aP,aU;uniform vec3 uCenter,uRight,uUp,uDepth;uniform vec2 uPan;uniform float uScale,uAspect,uRadius,uDeform;out vec3 world;out float value;void main(){vec3 a=aP+aU*uDeform-uCenter;world=a;value=length(aU);gl_Position=vec4(dot(a,uRight)/(uScale*uAspect)+uPan.x,dot(a,uUp)/uScale+uPan.y,-dot(a,uDepth)/(uRadius*4.),1.);}`));
  gl.attachShader(program,shader(gl.FRAGMENT_SHADER,`#version 300 es
precision highp float;in vec3 world;in float value;uniform float uMax;uniform bool uColor;out vec4 color;void main(){vec3 n=normalize(cross(dFdx(world),dFdy(world)));float shade=.80+.20*abs(dot(n,normalize(vec3(.3,-.5,.8))));float t=clamp(value/uMax,0.,1.);vec3 base=uColor?mix(vec3(.25,.1,.48),vec3(.96,.88,.18),t):vec3(.73,.78,.81);color=vec4(base*shade,1.);}`));
  gl.linkProgram(program);if(!gl.getProgramParameter(program,gl.LINK_STATUS))throw Error(gl.getProgramInfoLog(program));gl.enable(gl.DEPTH_TEST);gl.depthFunc(gl.LEQUAL);
  overlay=new FEBoundaryOverlay({canvas,gl,redraw:()=>draw(),initial:true});draw();
}
function frame(){
  const ca=Math.cos(az),sa=Math.sin(az),ce=Math.cos(el),se=Math.sin(el),cr=Math.cos(roll),sr=Math.sin(roll),r=[ca,sa,0],u=[-sa*se,ca*se,ce],right=r.map((v,i)=>v*cr+u[i]*sr),up=u.map((v,i)=>v*cr-r[i]*sr),aspect=canvas.width/canvas.height;
  return {id:caseId,sample:'peak',title:dataset.title,captionOffset:86,center:dataset._center,right,up,depth:[sa*ce,-ca*ce,se],pan,scale:dataset._radius*Math.max(.83,1.06/aspect)/zoom,radius:dataset._radius,deform:Number(document.getElementById('deform').value),viewport:{x:0,y:0,w:canvas.width,h:canvas.height}};
}
function draw(){
  if(!program||!overlay)return;
  if(!dataset){dataset=overlay.data(caseId);if(!dataset)return;overlay.mesh(dataset,dataset.samples.peak);const x=overlay.decode(dataset.xyz),u=overlay.decode(dataset.samples.peak.u);let lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity],maxU=0;
    for(let i=0;i<x.length;i++){lo[i%3]=Math.min(lo[i%3],x[i]);hi[i%3]=Math.max(hi[i%3],x[i]);}
    for(let i=0;i<u.length;i+=3)maxU=Math.max(maxU,Math.hypot(u[i],u[i+1],u[i+2]));dataset._maxU=maxU;dataset._center=lo.map((v,i)=>(v+hi[i])/2);dataset._radius=Math.hypot(...hi.map((v,i)=>v-lo[i]))/2;roll=dataset.kind==='axisym'?-1.15:0;
    const tri=overlay.decode(dataset.base_tri,Uint32Array);bodyTri=overlay.buffer(tri,gl.ELEMENT_ARRAY_BUFFER);bodyCount=tri.length;document.getElementById('title').textContent=dataset.title;document.getElementById('scope').textContent=dataset.kind==='axisym'?`原计算为轴对称二维 Q2：${dataset.original_nodes} 节点 / ${dataset.original_cells} 单元；绕轴网格仅供三维显示。`:'10 节点 / 1 个任意曲边四面体；保存的线性约束与构造节点力。';
  }
  const ratio=Math.min(devicePixelRatio||1,2);canvas.width=Math.round(canvas.clientWidth*ratio);canvas.height=Math.round(canvas.clientHeight*ratio);if(!canvas.width||!canvas.height)return;
  gl.useProgram(program);gl.viewport(0,0,canvas.width,canvas.height);gl.disable(gl.SCISSOR_TEST);gl.clearColor(.96,.978,.986,1);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);const f=frame();
  for(const k of ['center','right','up','depth'])gl.uniform3fv(gl.getUniformLocation(program,'u'+k[0].toUpperCase()+k.slice(1)),f[k]);gl.uniform2fv(gl.getUniformLocation(program,'uPan'),pan);
  for(const [key,value]of [['Scale',f.scale],['Aspect',canvas.width/canvas.height],['Radius',f.radius],['Deform',f.deform],['Max',Math.max(dataset._maxU,1e-20)]])gl.uniform1f(gl.getUniformLocation(program,'u'+key),value);
  gl.uniform1i(gl.getUniformLocation(program,'uColor'),Number(document.getElementById('color').value));
  for(const [name,b]of [['aP',dataset._pos],['aU',dataset.samples.peak._u]]){const loc=gl.getAttribLocation(program,name);gl.bindBuffer(gl.ARRAY_BUFFER,b);gl.enableVertexAttribArray(loc);gl.vertexAttribPointer(loc,3,gl.FLOAT,false,0,0);}
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,bodyTri);gl.drawElements(gl.TRIANGLES,bodyCount,gl.UNSIGNED_INT,0);overlay.render([f],program);
}
function reset(){az=.65;el=.35;zoom=1;pan=[0,0];draw();}
canvas.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.max(.3,Math.min(30,zoom*Math.exp(-e.deltaY*.0015)));draw();},{passive:false});
canvas.addEventListener('pointerdown',e=>{drag=[e.clientX,e.clientY,e.shiftKey||e.button===1];canvas.setPointerCapture(e.pointerId);});
canvas.addEventListener('pointermove',e=>{if(!drag)return;const dx=e.clientX-drag[0],dy=e.clientY-drag[1],shift=drag[2];drag=[e.clientX,e.clientY,shift];if(shift){pan[0]+=2*dx/canvas.clientWidth;pan[1]-=2*dy/canvas.clientHeight;}else{az+=dx*.006;el=Math.max(-1.5,Math.min(1.5,el+dy*.006));}draw();});
canvas.addEventListener('pointerup',()=>drag=null);canvas.addEventListener('pointercancel',()=>drag=null);canvas.addEventListener('dblclick',reset);
document.getElementById('reset').onclick=reset;document.getElementById('color').onchange=()=>draw();document.getElementById('deform').onchange=()=>draw();document.getElementById('case').onchange=e=>location.href='bc-inspector.html?case='+encodeURIComponent(e.target.value);new ResizeObserver(()=>draw()).observe(canvas);
try{init();}catch(error){document.getElementById('error').textContent=String(error);}
