/* Native surface fields, fixed external loads, and a declared shared display scale. */
(async()=>{'use strict';
const canvas=document.getElementById('scene'),gl=canvas.getContext('webgl2',{antialias:true,alpha:false,preserveDrawingBuffer:true});
async function load(url){const r=await fetch(url);if(!r.ok)throw Error(url);return url.endsWith('.gz')?new Response(r.body.pipeThrough(new DecompressionStream('gzip'))).json():r.json()}
function bytes(s,T){return new T(Uint8Array.from(atob(s),c=>c.charCodeAt(0)).buffer)}
function shader(type,src){const s=gl.createShader(type);gl.shaderSource(s,src);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(s));return s}
function normals(x,t){const n=new Float32Array(x.length);for(let k=0;k<t.length;k+=3){let a=t[k]*3,b=t[k+1]*3,c=t[k+2]*3;const u=[x[b]-x[a],x[b+1]-x[a+1],x[b+2]-x[a+2]],v=[x[c]-x[a],x[c+1]-x[a+1],x[c+2]-x[a+2]],p=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];for(const i of [a,b,c])for(let j=0;j<3;j++)n[i+j]+=p[j]}for(let i=0;i<n.length;i+=3){const d=Math.hypot(n[i],n[i+1],n[i+2])||1;for(let j=0;j<3;j++)n[i+j]/=d}return n}
try{
if(!gl)throw Error('WebGL2 unavailable');
const summary=await load('../assets/m01-checkpoint-v1/summary.json'),id=new URLSearchParams(location.search).get('id'),job=summary.renders.find(j=>j.id===id);if(!job)throw Error('Unknown render job');
const base='../assets/rod-review-visuals-v1/',manifest=await load(base+'manifest.json'),c=manifest.cases[job.case_index];
const [mesh,field]=await Promise.all([load(base+c.mesh),load(base+(job.source==='truth'?c.truth:`${job.source}-20-${String(job.case_index).padStart(2,'0')}.json.gz`))]);
const xyz=bytes(mesh.xyz,Float32Array),tri=bytes(mesh.tri,Uint32Array),u=bytes(field.u,Float32Array),vm=bytes(field.vm,Float32Array),pos=xyz.slice();for(let i=0;i<pos.length;i++)pos[i]+=u[i];
const palette=manifest.palettes.viridis,colors=new Float32Array(pos.length);for(let i=0;i<vm.length;i++){const value=job.field==='VM'?vm[i]:Math.hypot(u[3*i],u[3*i+1],u[3*i+2]),v=Math.max(0,Math.min(.999999,(value-job.range[0])/(job.range[1]-job.range[0])))*(palette.length-1),a=Math.floor(v),f=v-a;for(let j=0;j<3;j++)colors[i*3+j]=palette[a][j]*(1-f)+palette[a+1][j]*f}
const program=gl.createProgram();gl.attachShader(program,shader(gl.VERTEX_SHADER,`#version 300 es
precision highp float;in vec3 p,n,c;uniform vec3 center,right,up,depth;uniform float scale,aspect,ds;out vec3 N,C;void main(){vec3 x=p-center;gl_Position=vec4(dot(x,right)/(scale*aspect),dot(x,up)/scale,-dot(x,depth)/ds,1.);N=vec3(dot(n,right),dot(n,up),dot(n,depth));C=c;}`));
gl.attachShader(program,shader(gl.FRAGMENT_SHADER,`#version 300 es
precision highp float;in vec3 N,C;out vec4 color;void main(){vec3 n=normalize(N);if(n.z<0.)n=-n;vec3 light=normalize(vec3(-.35,.55,1.));float d=max(dot(n,light),0.);float s=pow(max(dot(n,normalize(light+vec3(0.,0.,1.))),0.),48.);color=vec4(C*(.68+.31*d)+vec3(.10*s),1.);}`));gl.linkProgram(program);if(!gl.getProgramParameter(program,gl.LINK_STATUS))throw Error(gl.getProgramInfoLog(program));gl.useProgram(program);
function attr(name,data){const b=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,b);gl.bufferData(gl.ARRAY_BUFFER,data,gl.STATIC_DRAW);const a=gl.getAttribLocation(program,name);gl.enableVertexAttribArray(a);gl.vertexAttribPointer(a,3,gl.FLOAT,false,0,0)}attr('p',pos);attr('n',normals(pos,tri));attr('c',colors);const indices=gl.createBuffer();gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,indices);gl.bufferData(gl.ELEMENT_ARRAY_BUFFER,tri,gl.STATIC_DRAW);
canvas.width=innerWidth;canvas.height=innerHeight;const width=canvas.width,height=canvas.height,aspect=width/height,az=-.98,el=1.05,ca=Math.cos(az),sa=Math.sin(az),ce=Math.cos(el),se=Math.sin(el),right=[ca,sa,0],up=[-sa*se,ca*se,ce],depth=[sa*ce,-ca*ce,se],center=job.bbox[0].map((v,i)=>(v+job.bbox[1][i])/2),corners=[];
for(const x of [0,1])for(const y of [0,1])for(const z of [0,1])corners.push([job.bbox[x][0]-center[0],job.bbox[y][1]-center[1],job.bbox[z][2]-center[2]]);
const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0),rx=Math.max(...corners.map(p=>Math.abs(dot(p,right)))),ry=Math.max(...corners.map(p=>Math.abs(dot(p,up)))),scale=Math.max(ry*1.18,rx/aspect*1.18),radius=Math.hypot(...job.bbox[1].map((v,i)=>v-job.bbox[0][i]))/2;
for(const [k,v] of Object.entries({center,right,up,depth}))gl.uniform3fv(gl.getUniformLocation(program,k),v);for(const [k,v] of Object.entries({scale,aspect,ds:scale*8}))gl.uniform1f(gl.getUniformLocation(program,k),v);
gl.viewport(0,0,width,height);gl.enable(gl.DEPTH_TEST);gl.depthFunc(gl.LEQUAL);gl.clearColor(1,1,1,1);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);gl.drawElements(gl.TRIANGLES,tri.length,gl.UNSIGNED_INT,0);
const painter=new RodAppliedLoads(gl),loads=await painter.load(c.id);painter.draw({center,right,up,depth,pan:[0,0],scale,radius,depthScale:scale*8,viewport:{x:0,y:0,w:width,h:height}},loads,job.source==='truth'?'FE':job.source+'_20');gl.finish();
canvas.setAttribute('aria-label','M01 rendered '+job.id);document.getElementById('status').hidden=true;document.body.dataset.ready=job.id;
}catch(e){document.getElementById('status').textContent='Render error: '+e.message;document.body.dataset.error=e.message;console.error(e)}
})();
