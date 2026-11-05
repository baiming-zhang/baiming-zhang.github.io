/* Sort only the M01 fixed-budget table; missing values always remain last. */
(()=>{'use strict';
 const table=document.getElementById('m01-budget-table');if(!table)return;
 const tbody=table.tBodies[0],status=document.getElementById('m01-budget-sort-status');
 const rows=[...tbody.querySelectorAll('tr[data-model]')];
 const state={step:6000,metric:'vm',direction:1,column:null};
 const names=['参数比 ≤1.05×','参数比 >1.05×','不适用或参数未确定'];
 const score=(row,metric)=>{const value=row.getAttribute('data-'+metric+state.step);return value===null?null:Number(value)};
 function compare(a,b){
   if(state.column!==null){const column=state.column;let av=a.cells[column].textContent.trim(),bv=b.cells[column].textContent.trim();const am=/^[-—–]?$/.test(av),bm=/^[-—–]?$/.test(bv);if(am||bm)return am===bm?Number(a.dataset.order)-Number(b.dataset.order):(am?1:-1);if(column===1||column===2){av=Number(av.replace(/[,×]/g,''));bv=Number(bv.replace(/[,×]/g,''));}return state.direction*(typeof av==='number'?av-bv:av.localeCompare(bv,'zh-CN',{numeric:true}))||Number(a.dataset.order)-Number(b.dataset.order);}
   const av=score(a,state.metric),bv=score(b,state.metric);
   if(av===null||bv===null)return av===bv?Number(a.dataset.order)-Number(b.dataset.order):(av===null?1:-1);
   const alternate=state.metric==='u'?'vm':'u';
   return state.direction*(av-bv)||state.direction*(score(a,alternate)-score(b,alternate))||Number(a.dataset.order)-Number(b.dataset.order);
 }
 function render(){
   const fragment=document.createDocumentFragment();
   names.forEach((name,group)=>{
     const header=document.createElement('tr');header.className='final-group';
     const cell=document.createElement('th');cell.colSpan=7;cell.scope='rowgroup';cell.textContent=name;header.append(cell);fragment.append(header);
     rows.filter(row=>Number(row.dataset.group)===group).sort(compare).forEach(row=>fragment.append(row));
   });tbody.replaceChildren(fragment);
   table.querySelectorAll('th[data-sort-step]').forEach(th=>{
     const selected=state.column===null&&Number(th.dataset.sortStep)===state.step;
     th.setAttribute('aria-sort',selected?(state.direction===1?'ascending':'descending'):'none');
     th.querySelectorAll('button').forEach(button=>{
       const active=selected&&(!button.dataset.metric||button.dataset.metric===state.metric);
       button.setAttribute('aria-pressed',String(active));
       const arrow=button.querySelector('.sort-arrow');if(arrow)arrow.textContent=active?(state.direction===1?' ↑':' ↓'):'';
     });
   });
   table.querySelectorAll('th[data-sort-column]').forEach(th=>{const active=Number(th.dataset.sortColumn)===state.column;th.setAttribute('aria-sort',active?(state.direction===1?'ascending':'descending'):'none');th.querySelector('button').setAttribute('aria-pressed',String(active));th.querySelector('.sort-arrow').textContent=active?(state.direction===1?' ↑':' ↓'):'';});
   status.textContent='当前排序：'+(state.column===null?state.step+'步 · '+state.metric.toUpperCase()+'误差':['配置','参数量','相对TS3参数','状态'][state.column])+(state.direction===1?'升序（低→高）':'降序（高→低）')+'；各参数组分别排序，缺失结果始终放在组末。';
 }
 table.querySelectorAll('button[data-step]').forEach(button=>button.addEventListener('click',()=>{
   const step=Number(button.dataset.step),metric=button.dataset.metric||state.metric;
   state.direction=state.column===null&&state.step===step&&state.metric===metric?-state.direction:1;
   state.column=null;state.step=step;state.metric=metric;render();
 }));
 [...table.tHead.rows[0].cells].slice(0,4).forEach((th,column)=>{const text=th.textContent;th.dataset.sortColumn=column;const button=document.createElement('button');button.type='button';button.className='sort-control';button.setAttribute('aria-label','按'+text+'排序');button.append(document.createTextNode(text));const arrow=document.createElement('span');arrow.className='sort-arrow';arrow.setAttribute('aria-hidden','true');button.append(arrow);th.textContent='';th.append(button);button.onclick=()=>{state.direction=state.column===column?-state.direction:1;state.column=column;render();};});render();
})();

/* Every other M01 table: click a column header; paired U/VM headers expose both metrics. */
(()=>{'use strict';
 document.querySelectorAll('#m01-final-results table,main.delivery-page table').forEach(table=>{
   if(table.id==='m01-budget-table')return;
   const heading=table.tHead?.rows[table.tHead.rows.length-1]||[...table.rows].find(row=>row.querySelector('th')&&!row.querySelector('td'));
   if(!heading)return;
   const columns=[...heading.cells];let selected=-1,metric=0,direction=1;
   const value=(cell,index)=>{let text=(cell?.textContent||'').trim();if(/^(?:-|—|–|N\/A|无|未完成|未开始)?$/.test(text))return null;if(index!==null&&index>=0)text=text.split(/[\/／]/)[index]||text;const number=text.replace(/,/g,'').match(/[-+]?\d*\.?\d+(?:e[-+]?\d+)?/i);return number&&/^[+\-]?(?:\d|\.)/.test(text)?Number(number[0]):text;};
   function sort(column,part){direction=selected===column&&metric===part?-direction:1;selected=column;metric=part;
     for(const body of table.tBodies){const rows=[...body.rows].filter(row=>row!==heading&&row.querySelector('td')&&row.cells.length===columns.length);rows.forEach((row,i)=>row.dataset.sortOriginal??=String(i));rows.sort((a,b)=>{const av=value(a.cells[column],part),bv=value(b.cells[column],part);if(av===null||bv===null)return av===bv?Number(a.dataset.sortOriginal)-Number(b.dataset.sortOriginal):(av===null?1:-1);return direction*(typeof av==='number'&&typeof bv==='number'?av-bv:String(av).localeCompare(String(bv),'zh-CN',{numeric:true}))||Number(a.dataset.sortOriginal)-Number(b.dataset.sortOriginal)});rows.forEach(row=>body.append(row));}
     columns.forEach((th,i)=>{th.setAttribute('aria-sort',i===column?(direction===1?'ascending':'descending'):'none');th.querySelectorAll('button').forEach(button=>{const active=i===column&&Number(button.dataset.part)===part;button.setAttribute('aria-pressed',String(active));button.querySelector('.sort-arrow').textContent=active?(direction===1?' ↑':' ↓'):'';});});
   }
   columns.forEach((th,column)=>{if(th.colSpan!==1)return;const label=th.textContent.trim();const paired=/U\s*[\/／]\s*(?:节点)?VM/i.test(label);th.textContent='';const make=(text,part)=>{const button=document.createElement('button');button.type='button';button.className='m01-any-sort';button.dataset.part=part;button.setAttribute('aria-label','按'+label+(paired?' '+text:'')+'排序');button.setAttribute('aria-pressed','false');button.append(document.createTextNode(text));const arrow=document.createElement('span');arrow.className='sort-arrow';arrow.setAttribute('aria-hidden','true');button.append(arrow);button.onclick=()=>sort(column,part);th.append(button);};if(paired){const prefix=label.replace(/U\s*[\/／]\s*(?:节点)?VM/i,'').trim();if(prefix){th.append(document.createTextNode(prefix));th.append(document.createElement('br'));}make('U',0);th.append(' / ');make('VM',1);}else make(label,-1);});
 });
})();
