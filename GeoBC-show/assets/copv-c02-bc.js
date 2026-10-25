'use strict';
// All anchors, masks and faces are exported from the executed C02 mesh/solve.
(() => {
  const bcButton=document.getElementById('bc-toggle');
  const loadButton=document.getElementById('load-toggle');
  const svg=document.getElementById('bc-annotations');
  const panel=document.getElementById('bc-panel');
  const content=document.getElementById('bc-content');
  let metadata=null, showBC=false, showLoads=false, loadMeshes=[], loadPromise=null, spatial=null;
  const baseDraw=draw;
  const xyzNames=['Ux','Uy','Uz'];
  const number=(v,n=3)=>Number(v).toFixed(n);
  const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
  const add=(a,b,k=1)=>a.map((v,i)=>v+k*b[i]);
  const norm=a=>Math.hypot(...a);
  const escape=s=>String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');

  function project(position, displacement=[0,0,0]) {
    const b=basis(),width=canvas.clientWidth,height=canvas.clientHeight,aspect=width/height;
    const extent=Math.max(b.u.reduce((s,v,i)=>s+Math.abs(v)*[125,125,320][i],0),b.r.reduce((s,v,i)=>s+Math.abs(v)*[125,125,320][i],0)/aspect)*1.18/zoom;
    const world=add(position,displacement,pressure*deformation);
    return {x:(dot(world,b.r)/(extent*aspect)+pan[0]+1)*width/2,y:(1-dot(world,b.u)/extent-pan[1])*height/2};
  }

  function renderSurfaces() {
    if(!showLoads || !loadMeshes.length)return;
    gl.useProgram(program);gl.disable(gl.DEPTH_TEST);gl.depthMask(false);
    gl.enable(gl.BLEND);gl.blendFunc(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA);
    for(const item of loadMeshes) {
      for(const [name,buff,size] of [['aPosition',item.p,3],['aU',item.u,3],['aStress',item.s,2]]) {
        gl.bindBuffer(gl.ARRAY_BUFFER,buff);const index=gl.getAttribLocation(program,name);
        gl.enableVertexAttribArray(index);gl.vertexAttribPointer(index,size,gl.FLOAT,false,0,0);
      }
      const color=item.color.slice();if(pressure===0)color[3]*=.35;
      gl.uniform4fv(gl.getUniformLocation(program,'uOverlay'),color);
      gl.uniform1i(gl.getUniformLocation(program,'uWire'),0);
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,item.tri);gl.drawElements(gl.TRIANGLES,item.count,gl.UNSIGNED_INT,0);
    }
    gl.uniform4f(gl.getUniformLocation(program,'uOverlay'),0,0,0,0);
    gl.disable(gl.BLEND);gl.depthMask(true);gl.enable(gl.DEPTH_TEST);
  }

  function renderGlyphs(){
    if(!spatial)spatial=new FESpatialGlyphs(gl);
    const b=basis(),aspect=canvas.width/canvas.height;
    const scale=Math.max(b.u.reduce((s,v,i)=>s+Math.abs(v)*[125,125,320][i],0),b.r.reduce((s,v,i)=>s+Math.abs(v)*[125,125,320][i],0)/aspect)*1.18/zoom;
    const frame={center:[0,0,0],right:b.r,up:b.u,depth:b.d,pan,scale,radius:375,depthScale:1500,deform:pressure*deformation,viewport:{x:0,y:0,w:canvas.width,h:canvas.height}},items=[];
    if(showLoads&&pressure>0){
      for(const sample of metadata.pressure.arrows){
        if(cut&&sample.position_mm[0]>.0001&&sample.position_mm[1]<-.0001)continue;
        items.push({p:sample.position_mm,u:sample.U_peak_mm,d:sample.direction,length:sample.length_mm,shift:-1,color:[1.,.36,.035]});
      }
      for(const cap of metadata.caps)items.push({p:cap.position_mm,u:cap.U_peak_mm,d:cap.direction,length:42,color:[.93,.18,.48]});
    }
    if(showBC)for(const point of metadata.constraints){
      items.push({p:point.position_mm,u:point.U_peak_mm,d:[0,0,1],length:2.6,shape:'ball',color:[.12,.65,1.]});
      for(const axis of point.locked_axes)items.push({p:point.position_mm,u:point.U_peak_mm,d:[0,1,2].map(i=>i===axis?1:0),length:18,shape:'support',color:[.12,.65,1.]});
    }
    spatial.draw(frame,items);
    const dpr=canvas.width/canvas.clientWidth;
    spatial.draw({...frame,center:[0,0,0],pan:[0,0],scale:1.5,radius:1,depthScale:4,deform:0,viewport:{x:10*dpr,y:80*dpr,w:90*dpr,h:90*dpr}},[0,1,2].map(i=>({p:[0,0,0],d:[0,1,2].map(j=>j===i?1:0),length:.9,color:[[.75,.19,.18],[.13,.52,.25],[.16,.38,.79]][i]})));
  }

  function renderAnnotations() {
    const visible=showBC||showLoads;svg.toggleAttribute('hidden',!visible);panel.hidden=!visible;
    canvas.parentElement.classList.toggle('bc-active',visible);
    if(!visible){svg.innerHTML='';return;}
    const width=canvas.clientWidth,height=canvas.clientHeight;
    svg.setAttribute('viewBox',`0 0 ${width} ${height}`);
    renderGlyphs();
    let markup='';
    if(showLoads && pressure>0) {
      for(const cap of metadata.caps) {
        const start=project(cap.position_mm,cap.U_peak_mm),end=project(add(cap.position_mm,cap.direction,42),cap.U_peak_mm);
        const value=cap.force_peak_N*pressure/metadata.peak_pressure_MPa;
        const x=Math.min(width-110,Math.max(10,end.x-35)),y=Math.max(18,end.y-12);
        markup+=`<text x="${x}" y="${y}" class="cap-label">${cap.id==='upper'?'上端':'下端'} ${number(value,2)} N</text>`;
      }
    }
    if(showBC) {
      const centerY=Math.max(65,Math.min(height-145,metadata.constraints.reduce((sum,p)=>sum+project(p.position_mm,p.U_peak_mm).y,0)/3));
      metadata.constraints.forEach((point,k)=>{
        const p=project(point.position_mm,point.U_peak_mm);
        if(p.x<-25||p.x>width+25||p.y<-25||p.y>height+25)return;
        const text=point.id+' · '+point.locked_axes.map(i=>xyzNames[i]).join('/')+'=0';
        const boxWidth=point.locked_axes.length===3?150:point.locked_axes.length===2?132:108;
        const x=Math.min(width-boxWidth-8,Math.max(8,p.x+(p.x<width*.55?40:-boxWidth-34)));
        const y=centerY+[-34,0,34][k];
        const lineEnd=x>p.x?x:x+boxWidth;
        const tooltip=`${point.id}: 原始节点 ${point.gmsh_node_id}; 坐标 [${point.position_mm.map(v=>number(v,4)).join(', ')}] mm; ${text}; 当前反力 [${point.reaction_peak_N.map(v=>(v*pressure).toExponential(3)).join(', ')}] N`;
        markup+=`<g><title>${escape(tooltip)}</title><path d="M${p.x} ${p.y} L${lineEnd} ${y}" fill="none" stroke="#2163b3" stroke-width="1.3"/><text x="${x+7}" y="${y+4}" class="constraint-label">${text}</text></g>`;
      });
    }
    const b=basis(),origin={x:55,y:height-125};
    markup+=`<text x="18" y="${origin.y-42}" class="axis-caption">全局坐标</text>`;
    for(const [i,color]of [[0,'#ba4545'],[1,'#32804b'],[2,'#3367ba']]){
      const x=origin.x+b.r[i]*27,y=origin.y-b.u[i]*27;
      markup+=`<text x="${x+4}" y="${y-4}" fill="${color}" class="axis-label">${'XYZ'[i]}</text>`;
    }
    svg.innerHTML=markup;
    let detail='';
    if(showBC) {
      detail+='<div class="bc-group"><b><i class="swatch constraint"></i>下端三点数值定位</b><div class="dof-list">';
      for(const point of metadata.constraints)detail+=`<span title="坐标 ${point.position_mm.map(v=>number(v,3)).join(', ')} mm；原始节点 ${point.gmsh_node_id}"><strong>${point.id}</strong> ${point.locked_axes.map(i=>xyzNames[i]).join(' / ')} = 0</span>`;
      const R=norm(metadata.gauge_reaction_sum_peak_N)*pressure/metadata.peak_pressure_MPa;
      detail+=`</div><p>全局位移分量；消除 6 个刚体模态。<br>墙面 / 铰链：无；结构界面理想粘结。</p><p>|ΣR| ≈ ${number(R,5)} N<br>对应已记录的端面积离散余量。</p></div>`;
    }
    if(showLoads) {
      const [upper,lower]=metadata.caps;
      detail+=`<div class="bc-group"><b><i class="swatch pressure"></i>内压 ${number(pressure,2)} MPa</b><p>橙色：内壁与端口内肩，${metadata.pressure.facet_count.toLocaleString()} 片。<br>t = −p n<sub>solid</sub>；外壁零表压。</p><b><i class="swatch cap"></i>闭口端等效合力</b><p>上端 +Z：${number(upper.force_peak_N*pressure,3)} N<br>下端 −Z：${number(lower.force_peak_N*pressure,3)} N<br>紫红色：实际环形端面。</p></div>`;
    }
    content.innerHTML=detail+'<p class="bc-footnote">显示开关不改变工况。三维实体标记透视显示，可配合剖面查看；箭头长度示意，数值随状态更新。</p>';
  }

  draw=function() {
    baseDraw();
    if(!metadata||!manifest||!program||!canvas.clientWidth||!canvas.clientHeight)return;
    renderSurfaces();renderAnnotations();
  };

  bcButton.onclick=()=>{
    showBC=!showBC;bcButton.classList.toggle('active',showBC);bcButton.setAttribute('aria-pressed',String(showBC));draw();
  };
  loadButton.onclick=async()=>{
    showLoads=!showLoads;loadButton.classList.toggle('active',showLoads);loadButton.setAttribute('aria-pressed',String(showLoads));
    if(showLoads&&!loadMeshes.length) {
      loadButton.textContent='载荷加载中…';
      try {
        if(!loadPromise)loadPromise=payload('bc-load-surfaces').then(data=>{loadMeshes=data.parts.map(part=>Object.assign(object(part),{color:part.color}));});
        await loadPromise;
      } catch(error) {
        loadPromise=null;showLoads=false;loadButton.classList.remove('active');loadButton.setAttribute('aria-pressed','false');document.getElementById('bc-summary').textContent='载荷位置载入失败，请重试。';
      } finally {loadButton.textContent='载荷';}
    }
    draw();
  };
  fetch('../assets/copv-c02/bc-loads.json').then(response=>{
    if(!response.ok)throw Error('BC metadata missing');return response.json();
  }).then(data=>{
    metadata=data;bcButton.disabled=false;loadButton.disabled=false;
    document.getElementById('bc-summary').textContent='约束：下端 3 点（3–2–1） · 载荷：内压 + 两端合力';draw();
  }).catch(()=>{document.getElementById('bc-summary').textContent='边界条件数据载入失败，请刷新。';});
})();
