'use strict';
const $=id=>document.getElementById(id);
const memory=new Map();
async function fetchData(file,binary=false){
  if(!memory.has(file)){
    const promise=fetch('./data/'+file).then(r=>{
      if(!r.ok)throw Error('HTTP '+r.status);
      return binary?r.arrayBuffer():r.json();
    }).catch(e=>{memory.delete(file);throw e});
    memory.set(file,promise);
  }
  return memory.get(file);
}
function option(select,value,text){
  const o=document.createElement('option');o.value=value;o.textContent=text;select.append(o);
}
function td(tr,value,tag='td'){
  const cell=document.createElement(tag);cell.textContent=value;tr.append(cell);return cell;
}
function fmt(value){
  if(value===null||value===undefined||value==='')return '—';
  if(typeof value!=='number')return String(value);
  if(!Number.isFinite(value))return '—';
  if(value===0)return '0';
  if(Math.abs(value)<.001)return value.toExponential(3);
  return value.toLocaleString('zh-CN',{maximumSignificantDigits:Number.isInteger(value)?15:7});
}
const views={
  arms:{name:'实验总表',start:2,ids:[0,1],note:'每组一次运行。正式误差仅基于该组有效题。'},
  pairs:{name:'配对变化与胜平负',start:1,ids:[0],note:'比较同模型的基础组与加约束组。误差胜平负仅统计双方均有效的题目。'},
  cases:{name:'逐题误差与计时',start:4,ids:[0,1,2,3],note:'每题的误差合并全部指定时次。计时范围可通过同一指标菜单查看。'},
  frames:{name:'逐时误差',start:6,ids:[0,1,2,4,5],note:'误差在当前题目的指定时刻计算。可切换正式、辅助和格式补算记录。'},
  physics:{name:'空间与物理诊断',start:6,ids:[0,1,2,4,5],note:'比较梯度、极值、积分和物理约束。不同题目适用的诊断不同。'},
  bounds:{name:'边界与界面',start:7,ids:[0,1,2,4,5,6],note:'每条记录对应一个壁面或界面。可查看对物理目标及离散参考的残差。'},
  aux:{name:'辅助运行与补算',start:5,ids:[0,1,2,3,4],note:'包括独立评分的细化运行和 Luna 基础组的格式敏感性检查。'},
  tele:{name:'回合与客户端事件',start:3,ids:[0,1,2],note:'精确内部推理时长无测量值。不同类型的客户端事件可以重叠。'},
  timers:{name:'原始计时记录',start:5,ids:[0,1,2,3,4,7,8],note:'可按阶段、文件名或原字段检索。条目含重复描述、时钟读数及先前快照。'},
  points:{name:'逐点偏差',start:5,ids:[0,1,2,3,4],note:'取一题的规定坐标，查看解析参考与选中模型的有符号误差。'}
};
const parts={'gradients':'梯度','all':'全部方向','edge_count':'边数','rmse':'RMSE','relative_l2':'相对 L2','maximum':'最大值','minimum':'最小值','prediction':'预测','reference':'参考','signed_amplitude_error':'振幅有符号误差','prediction_ties':'预测并列数','reference_ties':'参考并列数','first_prediction_xyz':'首个预测位置','first_reference_xyz':'首个参考位置','location_error_first_to_nearest':'首个位置至最近参考的距离','location_error_worst_to_nearest':'最大最近距离','integral':'积分','scope':'范围','reference_quadrature':'参考求积','signed_error':'有符号误差','absolute_error':'绝对误差','exact_global':'全域解析值','prediction_minus_global':'预测减全域解析值','quadrature_minus_global':'参考求积减全域解析值','physical_constraints':'物理约束','negative_count':'负值数','fraction_negative':'负值比例','max_negative_violation':'最大负向越界','above_one_count':'大于一的数量','max_above_one':'超过一的最大幅度','count':'样本数','target':'目标','prediction_mean':'预测均值','reference_mean':'参考均值','prediction_minus_target':'预测减物理目标','reference_minus_target':'参考减物理目标','prediction_minus_reference':'预测减参考','mean':'均值','max_abs':'最大绝对值'};
const label=s=>s.split('.').map(x=>parts[x]||x).join(' · ');
let page=0,version=0,shown=[],headers=[],models=[],arms;
const pageSize=15;
function summary(){
  const idx=Number($('measure').value),body=$('summary').tBodies[0];body.replaceChildren();
  models.forEach(m=>{
    const a=arms.rows.find(r=>r[0]===m&&r[1]==='基础'),b=arms.rows.find(r=>r[0]===m&&r[1]==='加约束');
    const tr=body.insertRow(),c=tr.insertCell(),btn=document.createElement('button');
    btn.type='button';btn.className='modelbutton';btn.textContent=m;btn.onclick=()=>jump(m,idx);c.append(btn);
    td(tr,fmt(a[idx]));td(tr,fmt(b[idx]));
    let ratio='—',cls='null';
    if(typeof a[idx]==='number'&&typeof b[idx]==='number'&&a[idx]!==0&&!(idx===18&&a[idx]<0)){
      const d=b[idx]/a[idx]-1;ratio=(d>0?'+':'')+(100*d).toLocaleString('zh-CN',{maximumFractionDigits:2})+'%';
      cls=d<0?'less':(d>0?'more':'');
    }
    td(tr,ratio).className='delta '+cls;td(tr,a[2]+'/10 · '+b[2]+'/10');
  });
  $('scale').textContent=(idx>=14?'误差均为无量纲数值。':'计时单位为秒，token 为记录数。')+'  变化率以基础提示为分母。零基数、负基数的偏差或缺测显示为横线。';
}
function jump(model,index){
  $('model').value=model;$('arm').value='';$('case').value='';$('kind').value='正式';$('search').value='';
  const mapping={11:13,12:14,13:15,14:6,16:7,17:8,18:9,19:10,20:11,21:12};
  const view=index in mapping?'cases':'arms';
  $('dataset').value=view;loadDetail(index in mapping?mapping[index]:index);
  $('detail').scrollIntoView({behavior:'smooth',block:'start'});
}
function renderPage(){
  const table=$('detailtable'),head=table.tHead,body=table.tBodies[0];
  head.replaceChildren();body.replaceChildren();
  const hr=head.insertRow();headers.forEach(h=>td(hr,label(h),'th'));
  const n=shown.length,last=Math.max(0,Math.ceil(n/pageSize)-1);page=Math.min(page,last);
  for(const r of shown.slice(page*pageSize,(page+1)*pageSize)){const tr=body.insertRow();r.forEach(v=>td(tr,fmt(v)))}
  if(!n){const c=td(body.insertRow(),'当前筛选没有适用记录。');c.colSpan=Math.max(1,headers.length)}
  $('count').textContent=n?'共 '+n.toLocaleString('zh-CN')+' 条 · 第 '+(page+1)+' / '+(last+1)+' 页':'0 条记录';
  $('previous').disabled=page===0;$('next').disabled=page===last;
}
async function loadDetail(wanted){
  const seq=++version,key=$('dataset').value,v=views[key],old=$('field').value;
  const point=key==='points',caseOK=!['arms','tele','pairs'].includes(key),hasType=['frames','bounds','physics','aux'].includes(key);
  $('arm').parentElement.hidden=key==='pairs';
  $('casebox').hidden=!caseOK;$('kindbox').hidden=!hasType;$('searchbox').hidden=key!=='timers';
  if(point&&!$('case').value)$('case').value='01';
  $('detailnote').textContent=v.note;$('grid').setAttribute('aria-busy','true');$('count').textContent='正在读取数据…';
  try{
    const payload=await fetchData(point?'point-cols.json':(key==='aux'?'extra-runs':key)+'.json');
    if(seq!==version)return;
    const cols=point?payload:payload.cols;
    const options=point?[[5,'解析参考值'],[6,'预测减参考']]:cols.map((x,i)=>[i,label(x)]).filter(([i])=>i>=v.start&&!(key==='timers'&&[7,8].includes(i)));
    const chosen=wanted??($('field').dataset.view===key?Number(old):({arms:14,pairs:15,cases:6,frames:7,physics:11,bounds:12,aux:7,tele:8,timers:5,points:6}[key]));
    $('field').replaceChildren();for(const [i,x] of options)option($('field'),i,x);
    $('field').value=String(options.some(([i])=>i===chosen)?chosen:options[0][0]);$('field').dataset.view=key;
    const field=Number($('field').value),model=$('model').value,arm=$('arm').value,cas=$('case').value,kind=$('kind').value;
    if(point){
      const bytes=await fetchData('points-'+cas+'.bin',true);if(seq!==version)return;
      const raw=new Float64Array(bytes),n=cols.length;
      const selected=field===5?[]:cols.map((x,i)=>[x,i]).filter(([x,i])=>i>=6&&(!model||x.startsWith(model+' ·'))&&(!arm||x.endsWith(' · '+arm))).map(x=>x[1]);
      const indices=[...v.ids,5,...selected];headers=indices.map(i=>cols[i]);shown=[];
      for(let k=0;k<raw.length;k+=n)shown.push(indices.map(i=>raw[k+i]));
    }else{
      const mi=key==='tele'?1:0,ai=key==='tele'?2:1,text=$('search').value.trim().toLowerCase();
      const rows=payload.rows.filter(r=>(!model||r[mi]===model)&&(key==='pairs'||!arm||r[ai]===arm)&&(!caseOK||!cas||r[2]===cas)&&(!hasType||!kind||r[3]===kind)&&(key!=='timers'||!text||r.some(x=>String(x??'').toLowerCase().includes(text))));
      const indices=key==='timers'?[0,1,2,3,field,4,7,8]:[...v.ids,field];headers=indices.map(i=>cols[i]);shown=rows.map(r=>indices.map(i=>r[i]));
    }
    page=0;renderPage();
  }catch(e){if(seq===version){shown=[];headers=[];renderPage();$('count').textContent='数据读取失败，请刷新重试。'}}
  finally{if(seq===version)$('grid').setAttribute('aria-busy','false')}
}
async function init(){
  arms=await fetchData('arms.json');models=[...new Set(arms.rows.map(r=>r[0]))];
  for(const i of [14,15,16,17,18,19,20,21,3,4,5,6,7,8,9,10,11,12,13,2,22,23,24])option($('measure'),i,arms.cols[i]);
  for(const [k,v] of Object.entries(views))option($('dataset'),k,v.name);
  option($('model'),'','全部模型');models.forEach(m=>option($('model'),m,m));
  option($('arm'),'','两组');['基础','加约束'].forEach(a=>option($('arm'),a,a));
  option($('case'),'','全部题目');
  const names=['周期平流衰减','一维热扩散','绝热矩形扩散','三维各向异性输运','三维受迫扩散','入口通量扩散','复合板稳态导热','黏性 Burgers','旋转扩散','三维点源扩散'];
  names.forEach((name,i)=>{const k=String(i+1).padStart(2,'0');option($('case'),k,k+' · '+name)});
  option($('kind'),'','全部结果');['正式','辅助','格式补算'].forEach(a=>option($('kind'),a,a));
  $('dataset').value='cases';$('model').value=models[0];$('kind').value='正式';
  $('measure').onchange=summary;
  $('dataset').onchange=()=>{if($('dataset').value==='aux')$('kind').value='';else if(['frames','physics','bounds'].includes($('dataset').value))$('kind').value='正式';loadDetail()};
  for(const k of ['model','arm','case','field','kind'])$(k).onchange=()=>loadDetail();
  $('search').oninput=()=>loadDetail();
  $('previous').onclick=()=>{page--;renderPage()};$('next').onclick=()=>{page++;renderPage()};
  summary();await loadDetail();
}
init().catch(()=>{$('scale').textContent='无法读取总表。请刷新页面，或下载完整 Excel。'});
