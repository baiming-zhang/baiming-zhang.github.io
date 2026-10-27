/* Same closed 3D arrow geometry and pressure/body colors as C04. */
window.RodAppliedLoads=class {
  constructor(gl){this.gl=gl;this.cache=new Map();this.spatial=new FESpatialGlyphs(gl);this.manifest=fetch('../assets/rod-load-glyphs-v1/manifest.json').then(r=>{if(!r.ok)throw Error('载荷清单未载入');return r.json()})}
  async load(id){
    if(!this.cache.has(id))this.cache.set(id,this.manifest.then(async m=>{const item=m.cases[id];if(!item)throw Error('此工况没有核验过的载荷标记：'+id);const r=await fetch('../assets/rod-load-glyphs-v1/'+item.asset);if(!r.ok)throw Error('原始载荷位置未载入');const data=await new Response(r.body.pipeThrough(new DecompressionStream('gzip'))).json();if(data.sample_id!==id)throw Error('载荷样本身份不匹配');return data}));
    return this.cache.get(id)
  }
  draw(frame,data,variant='FE'){
    if(!data)return;const reference=variant==='reference',u=reference?null:variant==='FE'?data.U_FE_mm:data.U_predictions_mm[variant];if(!reference&&!u)throw Error('缺少当前预测场的载荷位置映射');
    const items=data.glyphs.map((g,i)=>({p:g.p,u:reference?[0,0,0]:u[i],d:g.d,shift:g.shift,length:frame.radius*(g.kind==='pressure'?.14:.24),color:g.kind==='body'?[.93,.18,.48]:[1.,.36,.035]}));
    this.spatial.draw({...frame,deform:reference?0:1},items)
  }
};
