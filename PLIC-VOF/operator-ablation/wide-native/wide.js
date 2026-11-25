
(async()=>{
const $=id=>document.getElementById(id),esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
try{
 const d=window.WIDE_DATA||await fetch('data/reader-data.json').then(r=>{if(!r.ok)throw Error('HTTP '+r.status);return r.json();});
 const a=d.accuracy.rows.map(r=>Object.fromEntries(d.accuracy.fields.map((f,i)=>[f,r[i]]))),t=d.timing,labels=d.labels;
 const number=v=>{if(v===''||v===null||v===undefined)return '—';const n=Number(v);if(!Number.isFinite(n))return '—';return Math.abs(n)>0&&Math.abs(n)<.0001?n.toExponential(2):Number(n.toPrecision(5)).toLocaleString('en-US',{maximumFractionDigits:8});};
 const pm=(m,s)=>m===''?'无完整输出':number(m)+(s===''?'':' ± '+number(s));
 const name=r=>{let m=r.method||r.group_id;if(m.startsWith('aphros_native'))return 'APHROS：'+(r.operator==='flux'?'PLIC通量':'粒子串曲率');if(m.startsWith('simliquid_native'))return 'SimLiquid：'+(r.operator==='flux'?'PLIC通量':'粒子串曲率');if(m.includes('original_frozen'))return '此前训练的3×3联合网络';if(m.includes('height_function'))return '独立高度函数参考';if(m==='youngs_plic_flux')return '独立Youngs参考';const q=m.match(/(control|enhanced)_(?:s)?(\d+)(?:x|$)/);if(q)return '重训'+q[2]+'×'+q[2]+' · '+labels[q[1]];return m;};
 function options(id,v,selected){$(id).innerHTML=v.map(x=>'<option value="'+esc(x)+'">'+esc(labels[x]||x)+'</option>').join('');$(id).value=selected;}
 options('layer',['CPU_FP64','CPU_FP32','GPU_FP32'],'CPU_FP64');options('geometry',['independent_noncubic','physical_graph'],'independent_noncubic');options('grid',['32','64','128','256','512'],'128');options('metric',['physical_kappa','q','face_average_flux_density','ratio','cell_transport_fraction'],'physical_kappa');options('stratum',['all','plane','weak','curved','high_order'],'all');options('support',['all_methods_common','all_methods_native_finite_common'],'all_methods_native_finite_common');
 const layer=r=>r.method.startsWith('CPU_FP64:')?'CPU_FP64':r.method.startsWith('CPU_FP32:')?'CPU_FP32':r.method.startsWith('GPU_FP32:')||r.method==='original_frozen'?'GPU_FP32':'reference';
 const key=r=>[r.fixture_id,r.operator,r.metric,r.stratum,r.N,r.CFL,r.method].join('|');const original=new Map(a.filter(r=>r.subset==='all_methods_common').map(r=>[key(r),r]));
 function drawAccuracy(){const rr=a.filter(r=>r.subset===$('support').value&&r.fixture_id===$('geometry').value&&r.N===$('grid').value&&r.metric===$('metric').value&&r.stratum===$('stratum').value&&(layer(r)===$('layer').value||layer(r)==='reference')&&r.samples!=='0');$('accuracy').innerHTML=rr.length?rr.map(r=>{const o=original.get(key(r))||r;return '<tr><td>'+esc(name(r))+'</td><td>'+pm(r.mean,r.std)+'</td><td>'+esc(r.samples)+'</td><td>'+number(100*Number(o.prediction_coverage_min))+'%</td><td>'+esc(o.samples)+'</td></tr>';}).join(''):'<tr><td colspan="5">该几何集没有这一分层。</td></tr>';}
 ['layer','geometry','grid','metric','stratum','support'].forEach(id=>$(id).addEventListener('change',drawAccuracy));drawAccuracy();
 options('time-geometry',['independent_noncubic','physical_graph'],'independent_noncubic');options('time-grid',['128','512'],'512');options('time-operator',['flux','curvature'],'flux');options('time-batch',['1','32','256','1024','4096'],'4096');
 function drawTiming(){const rr=t.filter(r=>r.fixture_id===$('time-geometry').value&&r.N===$('time-grid').value&&r.operator===$('time-operator').value&&r.batch_size===$('time-batch').value);$('timing').innerHTML=rr.map(r=>'<tr><td>'+esc(name(r))+'</td><td>'+pm(r.error_mean,r.error_std)+'</td><td>'+pm(r.wall_ms_mean,r.wall_ms_std)+'</td><td>'+number(r.p05_mean_ms)+'–'+number(r.p95_mean_ms)+'</td><td>'+number(100*Number(r.full_task_success_rate_min))+'%</td></tr>').join('');}
 ['time-geometry','time-grid','time-operator','time-batch'].forEach(id=>$(id).addEventListener('change',drawTiming));drawTiming();$('load-status').textContent='全部精度表与同CPU计时表已载入。';
}catch(e){$('load-status').textContent='表格暂未载入，请刷新，或下载下方完整数据包。';console.error(e);}
})();
