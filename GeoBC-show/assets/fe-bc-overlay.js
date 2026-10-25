'use strict';
// Reusable view-only overlay. Every position/U/normal/load is exported from
// the frozen inputs; the adapter supplies the original viewer's exact camera.
window.FEBoundaryOverlay=class {
  constructor(options) {
    Object.assign(this,options);this.cache=new Map();this.bc=false;this.loads=false;
    this.host=this.canvas.parentElement;this.svg=document.createElementNS('http://www.w3.org/2000/svg','svg');
    this.svg.classList.add('fe-bc-svg');this.svg.setAttribute('aria-label','实际边界条件与载荷');this.svg.hidden=true;this.host.append(this.svg);
    const toolbar=document.querySelector('.toolbar');
    for(const [key,label]of [['bc','边界条件'],['loads','载荷']]) {
      const button=document.createElement('button');button.type='button';button.textContent=label;button.setAttribute('aria-pressed','false');button.className='fe-bc-button';
      button.onclick=()=>{this[key]=!this[key];button.setAttribute('aria-pressed',String(this[key]));this.redraw();};toolbar.append(button);this[key+'Button']=button;
    }
    this.info=document.createElement('details');this.info.className='fe-bc-info';this.info.hidden=true;
    this.info.innerHTML='<summary>工况 / 约束说明</summary><div class="fe-bc-cards"></div><p>蓝色：约束；橙色：表面压力 / 节点力；紫红色：惯性体载荷（如有）。标记透视显示，箭头长度示意；显示开关不改变工况。</p>';
    toolbar.after(this.info);
    const params=new URLSearchParams(location.search);this.bc=params.get('bc')==='1'||!!options.initial;this.loads=params.get('loads')==='1'||!!options.initial;
    this.bcButton.setAttribute('aria-pressed',String(this.bc));this.loadsButton.setAttribute('aria-pressed',String(this.loads));
  }
  decode(s,Type=Float32Array){return new Type(Uint8Array.from(atob(s),c=>c.charCodeAt(0)).buffer);}
  buffer(data,target){const g=this.gl,b=g.createBuffer();g.bindBuffer(target||g.ARRAY_BUFFER,b);g.bufferData(target||g.ARRAY_BUFFER,data,g.STATIC_DRAW);return b;}
  data(id) {
    if(this.cache.has(id))return this.cache.get(id).data;
    const item={data:null};this.cache.set(id,item);
    fetch('../assets/fe-bc/'+id+'.json.gz').then(r=>{if(!r.ok)throw Error('位置数据加载失败');return new Response(r.body.pipeThrough(new DecompressionStream('gzip'))).json();}).then(data=>{item.data=data;this.redraw();}).catch(error=>{item.error=String(error);this.redraw();});return null;
  }
  shader(type,s){const g=this.gl,q=g.createShader(type);g.shaderSource(q,s);g.compileShader(q);if(!g.getShaderParameter(q,g.COMPILE_STATUS))throw Error(g.getShaderInfoLog(q));return q;}
  initGL() {
    if(this.overlayProgram)return;
    const g=this.gl,p=g.createProgram();g.attachShader(p,this.shader(g.VERTEX_SHADER,`#version 300 es
in vec3 aP,aU;uniform vec3 uCenter,uRight,uUp,uDepth;uniform vec2 uPan;uniform float uScale,uAspect,uRadius,uDeform;out vec3 ref;void main(){ref=aP;vec3 a=aP+aU*uDeform-uCenter;gl_Position=vec4(dot(a,uRight)/(uScale*uAspect)+uPan.x,dot(a,uUp)/uScale+uPan.y,-dot(a,uDepth)/(uRadius*4.),1.);}`));
    g.attachShader(p,this.shader(g.FRAGMENT_SHADER,`#version 300 es
precision highp float;in vec3 ref;uniform bool uMask;uniform vec2 uDirection;uniform vec4 uColor;out vec4 color;void main(){if(uMask&&dot(normalize(ref.xy),uDirection)<=0.)discard;color=uColor;}`));
    g.linkProgram(p);if(!g.getProgramParameter(p,g.LINK_STATUS))throw Error(g.getProgramInfoLog(p));this.overlayProgram=p;
  }
  mesh(data,sample) {
    if(!data._pos){data._pos=this.buffer(this.decode(data.xyz));data._groups=data.groups.map(part=>({...part,b:this.buffer(this.decode(part.tri,Uint32Array),this.gl.ELEMENT_ARRAY_BUFFER),count:this.decode(part.tri,Uint32Array).length}));}
    if(!sample._u)sample._u=this.buffer(this.decode(sample.u));
  }
  project(frame,p,u=[0,0,0]) {
    const a=p.map((x,i)=>x+(frame.deform||0)*u[i]-frame.center[i]);const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),v=frame.viewport;
    const sx=this.canvas.width/this.canvas.clientWidth,sy=this.canvas.height/this.canvas.clientHeight;
    return {x:(v.x+(dot(a,frame.right)/(frame.scale*(v.w/v.h))+frame.pan[0]+1)*v.w/2)/sx,y:(this.canvas.height-v.y-v.h+(1-dot(a,frame.up)/frame.scale-frame.pan[1])*v.h/2)/sy};
  }
  escape(s){return String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');}
  fmt(v){return v==null?'—':Number(v).toLocaleString('zh-CN',{maximumFractionDigits:3});}
  vector(a){return a?a.map(x=>this.fmt(x)).join(', '):'—';}
  render(frames,baseProgram) {
    const enabled=this.bc||this.loads;this.svg.toggleAttribute('hidden',!enabled);this.info.hidden=!enabled;
    if(!enabled){this.svg.innerHTML='';return;}
    if(!this.gl||!baseProgram||!this.canvas.clientHeight||!this.canvas.clientWidth)return;
    this.initGL();const g=this.gl,p=this.overlayProgram;let drawings='',cards=[];const unique=new Set();
    this.svg.setAttribute('viewBox',`0 0 ${this.canvas.clientWidth} ${this.canvas.clientHeight}`);
    let defs='<defs><marker id="fe-pressure-head" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7" markerHeight="7" orient="auto" markerUnits="userSpaceOnUse"><path d="M0 0L8 4L0 8Z" fill="#e05b18"/></marker><marker id="fe-body-head" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7" markerHeight="7" orient="auto" markerUnits="userSpaceOnUse"><path d="M0 0L8 4L0 8Z" fill="#bb285b"/></marker>';
    g.useProgram(p);g.disable(g.DEPTH_TEST);g.depthMask(false);g.enable(g.BLEND);g.blendFunc(g.SRC_ALPHA,g.ONE_MINUS_SRC_ALPHA);g.enable(g.SCISSOR_TEST);
    for(let fi=0;fi<frames.length;fi++) {
      const f=frames[fi],d=this.data(f.id);if(!d)continue;const s=d.samples[f.sample];if(!s)continue;
      this.mesh(d,s);const v=f.viewport;g.viewport(v.x,v.y,v.w,v.h);g.scissor(v.x,v.y,v.w,v.h);
      const u=(name,value)=>{const loc=g.getUniformLocation(p,name);if(Array.isArray(value)){if(value.length===3)g.uniform3fv(loc,value);else g.uniform2fv(loc,value);}else g.uniform1f(loc,value);};
      for(const name of ['center','right','up','depth','pan'])u('u'+name[0].toUpperCase()+name.slice(1),f[name]);
      u('uScale',f.scale);u('uAspect',v.w/v.h);u('uRadius',f.radius);u('uDeform',f.deform||0);u('uDirection',s.pressure_direction||[1,0]);
      for(const [name,buffer]of [['aP',d._pos],['aU',s._u]]){const loc=g.getAttribLocation(p,name);g.bindBuffer(g.ARRAY_BUFFER,buffer);g.enableVertexAttribArray(loc);g.vertexAttribPointer(loc,3,g.FLOAT,false,0,0);}
      for(const group of d._groups) {
        if(group.kind==='constraint'?!this.bc:!this.loads)continue;
        const color=group.kind==='constraint'?[.12,.37,.80,.24]:[1.,.34,.03,.23];
        g.uniform4fv(g.getUniformLocation(p,'uColor'),color);g.uniform1i(g.getUniformLocation(p,'uMask'),group.kind==='pressure'&&d.kind==='rod'?1:0);
        g.bindBuffer(g.ELEMENT_ARRAY_BUFFER,group.b);g.drawElements(g.TRIANGLES,group.count,g.UNSIGNED_INT,0);
      }
      const sx=this.canvas.width/this.canvas.clientWidth,sy=this.canvas.height/this.canvas.clientHeight,left=v.x/sx,top=(this.canvas.height-v.y-v.h)/sy,w=v.w/sx,h=v.h/sy;
      defs+=`<clipPath id="fe-bc-clip-${fi}"><rect x="${left}" y="${top}" width="${w}" height="${h}"/></clipPath>`;
      drawings+=`<g clip-path="url(#fe-bc-clip-${fi})">`;
      for(const item of s.glyphs) {
        if(item.kind==='constraint'?!this.bc:!this.loads)continue;
        const at=this.project(f,item.p,item.u),length=f.radius*(item.kind==='constraint'?.065:.13);
        if(item.kind==='constraint') {
          const end=this.project(f,item.p.map((x,i)=>x+item.d[i]*length),item.u),dx=end.x-at.x,dy=end.y-at.y,mag=Math.hypot(dx,dy),tx=mag>1?-dy/mag*3:3,ty=mag>1?dx/mag*3:0;
          drawings+=`<path d="M${at.x} ${at.y}L${end.x} ${end.y} M${end.x-tx} ${end.y-ty}L${end.x+tx} ${end.y+ty}" stroke="white" stroke-width="4" fill="none"/><path d="M${at.x} ${at.y}L${end.x} ${end.y} M${end.x-tx} ${end.y-ty}L${end.x+tx} ${end.y+ty}" stroke="#256ac0" stroke-width="1.5" fill="none"/><circle cx="${at.x}" cy="${at.y}" r="2.4" fill="#256ac0"/>`;
        } else {
          const start=item.kind==='pressure'?this.project(f,item.p.map((x,i)=>x-item.d[i]*length),item.u):at;
          const end=item.kind==='pressure'?at:this.project(f,item.p.map((x,i)=>x+item.d[i]*length),item.u);
          const color=item.kind==='body'?'#bb285b':'#e05b18',marker=item.kind==='body'?'fe-body-head':'fe-pressure-head';
          if(Math.hypot(end.x-start.x,end.y-start.y)<2)drawings+=`<circle cx="${at.x}" cy="${at.y}" r="2" fill="${color}"/>`;
          else drawings+=`<path d="M${start.x} ${start.y}L${end.x} ${end.y}" stroke="white" stroke-width="4"/><path d="M${start.x} ${start.y}L${end.x} ${end.y}" stroke="${color}" stroke-width="1.6" marker-end="url(#${marker})"/>`;
        }
      }
      let lines=[];
      if(this.bc)lines.push(d.kind==='rod'?'径向 uᵣ=0；〈uθ〉A=〈uz〉A=0':d.kind==='axisym'?'底端 uz=0；ur 自由，上端自由':'nᵢ·Uᵢ=0 + 两条均值线性方程');
      if(this.loads)lines.push(d.kind==='rod'?`销孔 p₀=${this.fmt(s.p0_MPa)} MPa；紫红色为惯性体载荷`:d.kind==='axisym'?'内压 1 MPa；无闭口端合力':'构造节点力（数值核对输入）');
      if(f.showDeform)lines.push(f.deform?'位移显示 ×'+Number(f.deform).toFixed(2):'参考几何（未放大位移）');
      lines.forEach((text,k)=>drawings+=`<text x="${left+10}" y="${top+(f.captionOffset||58)+15*k}" class="fe-bc-caption">${this.escape(text)}</text>`);
      drawings+='</g>';
      const key=f.id+'__'+f.sample;
      if(!unique.has(key)) {
        unique.add(key);let html=`<article><b>${this.escape(f.title||f.id)}</b> <span class="${s.state==='input_only'?'input-only':''}">${s.state==='input_only'?'规定满载输入 · 未接受解':'已接受满载工况'}</span>`;
        if(this.bc)html+='<p>'+d.boundary_text.map(x=>this.escape(x)).join('<br>')+`<br>约束节点：${d.support_nodes}</p>`;
        if(this.loads) {
          if(d.kind==='rod')html+=`<p>销孔半周径向压力：p₀=${this.fmt(s.p0_MPa)} MPa，κ=${this.fmt(s.skew)}<br>F销孔=[${this.vector(s.pin_force_N)}] N<br>F惯性体载荷=[${this.vector(s.body_force_N)}] N<br>发动机净气压=${this.fmt(s.snapshot?.net_pressure_MPa)} MPa（载荷参数，不是均布销孔压力）</p>`;
          else if(d.kind==='axisym')html+='<p>内表面 1 MPa；上端牵引自由。<br>开放通口接入试验，不添加 C02 闭口端力。</p>';
          else html+='<p>橙色为保存输入中的逐节点力，方向与数值由 force_N 读取；完整分量见位置数据文件。</p>';
        }
        cards.push(html+'</article>');
      }
    }
    g.disable(g.BLEND);g.depthMask(true);g.enable(g.DEPTH_TEST);g.disable(g.SCISSOR_TEST);g.useProgram(baseProgram);
    this.svg.innerHTML=defs+'</defs>'+drawings;
    this.info.querySelector('summary').textContent=`工况 / 约束说明 · ${this.bc?'蓝色：约束 ':''}${this.loads?'橙色：表面/节点载荷；紫红色：体载荷':''}（展开数值）`;
    this.info.querySelector('.fe-bc-cards').innerHTML=cards.join('')||'正在读取已保存的工况位置数据…';
  }
};
