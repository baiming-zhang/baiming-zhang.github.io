
(() => {
 const data=window.STUDY_DATA, $=id=>document.getElementById(id), labels=data.labels;
 const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
 const num=v=>{if(v===''||v===null||v===undefined)return '—';v=Number(v);if(!Number.isFinite(v))return '—';if(v===0)return '0';return Math.abs(v)<.0001?v.toExponential(2):Number(v.toPrecision(4)).toLocaleString('en-US',{maximumFractionDigits:7});};
 const pm=(m,s)=>s===''?num(m):num(m)+' ± '+num(s);
 function options(id,values,chosen){$(id).innerHTML=values.map(v=>'<option value="'+esc(v)+'">'+esc(labels[v]||v)+'</option>').join('');$(id).value=chosen;}
 options('a-stage',['12000','1500'],'12000');options('a-fixture',['independent_noncubic','physical_graph'],'independent_noncubic');options('a-N',['32','64','128','256','512'],'128');options('a-stratum',[...new Set(data.accuracy.map(r=>r.stratum))],'all');options('a-metric',['physical_kappa','q','face_average_flux_density','ratio','cell_transport_fraction'],'physical_kappa');
 function drawAccuracy(){const list=data.accuracy.filter(r=>r.stage===$('a-stage').value&&r.fixture_id===$('a-fixture').value&&r.N===$('a-N').value&&r.stratum===$('a-stratum').value&&r.metric===$('a-metric').value);$('accuracy-body').innerHTML=list.length?list.map(r=>'<tr><td>'+esc(labels[r.method]||r.method)+'</td><td>'+pm(r.mean,r.std)+'</td><td>'+esc(r.seeds||'单个固定参考')+'</td><td>'+esc(r.samples)+'</td><td>'+num(100*Number(r.coverage))+'%</td></tr>').join(''):'<tr><td colspan="5">该几何测试集没有这一分层，请更换分层。</td></tr>';$('accuracy-status').textContent=list.length+'种方法；三种子结果以均值 ± 样本标准差表示。';}
 ['a-stage','a-fixture','a-N','a-stratum','a-metric'].forEach(id=>$(id).addEventListener('change',drawAccuracy));drawAccuracy();
 options('t-stage',['12000','1500'],'12000');options('t-op',['flux','curvature'],'flux');options('t-batch',['1','32','256','1024','4096'],'4096');
 function drawTiming(){const list=data.timing.filter(r=>r.stage===$('t-stage').value&&r.operator===$('t-op').value&&r.batch===$('t-batch').value);$('timing-body').innerHTML=list.map(r=>'<tr><td>'+esc(labels[r.method]||r.method)+'</td><td>'+Number(r.mean).toFixed(1)+' ± '+Number(r.std).toFixed(1)+'</td><td>'+Number(r.p05).toFixed(1)+'–'+Number(r.p95).toFixed(1)+'</td><td>'+Number(r.per_stencil).toFixed(3)+'</td></tr>').join('');}
 ['t-stage','t-op','t-batch'].forEach(id=>$(id).addEventListener('change',drawTiming));drawTiming();
})();
