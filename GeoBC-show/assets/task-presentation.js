/* Model-first progressive presentation. Move original nodes; retain all content and handlers. */
(()=>{'use strict';
 const main=document.querySelector('main.delivery-page'),task=document.body.dataset.task;
 if(!main||!task||main.dataset.presentation==='model-first')return;
 const title=main.querySelector(':scope > .delivery-title');if(!title)return;
 const originalChildren=[...main.children],originalImages=[...main.querySelectorAll('img[src]')];
 const originalFrames=[...main.querySelectorAll('iframe[src]')];
 const discussions=[...main.querySelectorAll('section.discussion')];
 const looseSections=[...document.querySelectorAll('body > section')].filter(n=>!n.matches('.discussion'));
 const originalReview=main.querySelector('#public-review')?.textContent.trim()||'';
 const originalSummary=main.querySelector('#public-summary')?.textContent.trim()||'';
 const status=title.querySelector('.delivery-status')?.textContent.trim()||'';
 let detail=originalChildren.find(n=>n.matches('details.task-details'));
 if(!detail){detail=document.createElement('details');detail.className='task-details';const summary=document.createElement('summary');summary.textContent='展开详情：参数、数据、下载与完整审阅';detail.append(summary);main.append(detail);}
 detail.classList.add('task-presentation-details');detail.removeAttribute('open');
 let content=detail.querySelector(':scope > .details-content');
 if(!content){content=document.createElement('div');content.className='details-content';detail.append(content);}
 let stage=originalChildren.find(n=>n.matches('.delivery-stage')&&n.querySelector('iframe'));
 if(task==='M01')stage=main.querySelector('section[data-model-case="0095"]')||stage;
 // M01 reuses the already delivered native FE viewer, clearly labelled as FE context.
 if(!stage&&task==='M01'){
  stage=document.createElement('section');stage.className='delivery-stage task-model-context';
  const frame=document.createElement('iframe');frame.src='viewers/r05-models-20.html?v=loads-v1';frame.title='原始FE连杆场 · 已交付20例可旋转模型';frame.loading='eager';frame.allowFullscreen=true;stage.append(frame);
 }
 const frame=stage?.querySelector('iframe');let model=null;
 if(frame){
  model=document.createElement('section');model.className='task-priority-model';
  const tools=document.createElement('div');tools.className='task-model-tools';
  const label=document.createElement('strong');label.textContent=task==='M01'?(frame.getAttribute('src').includes('m01-case0095')?'0095 · 真值 / TS3 / MGN / 误差':'原始FE模型 · 交互查看'):'模型 · 旋转、缩放与场查看';tools.append(label);
  const choices=[{src:frame.getAttribute('src'),title:frame.title||'当前模型'}];
  originalFrames.forEach(f=>{const src=f.getAttribute('src');if(src&&!choices.some(x=>x.src===src))choices.push({src,title:f.title||'其他模型视图'});});
  if(task==='M01'&&!choices.some(x=>x.src.includes('r05-models-20')))choices.push({src:'viewers/r05-models-20.html?v=loads-v1',title:'原始FE20例模型 · 上下文参考'});
  if(task==='R07'){
   const slices=main.querySelector('a[href*="viewers/ts3-slices.html"]');
   if(slices)choices.push({src:slices.getAttribute('href'),title:'64 slices空间分配 · 可旋转模型'});
  }
  const fullModelURL=src=>{if(src.includes('m01-case0095')){const u=new URL(src,location.href);u.searchParams.delete('embed');return u.href;}return src;};
  const full=document.createElement('a');full.href=fullModelURL(frame.getAttribute('src'));full.target='_blank';full.rel='noopener';full.textContent='全屏打开 ↗';
  if(choices.length>1){const select=document.createElement('select');select.setAttribute('aria-label','选择模型视图');choices.forEach(v=>{const o=document.createElement('option');o.value=v.src;o.textContent=v.title;select.append(o);});select.addEventListener('change',()=>{frame.src=select.value;full.href=fullModelURL(select.value);frame.title=select.selectedOptions[0].textContent;if(task==='M01')label.textContent=select.selectedOptions[0].textContent;});tools.append(select);}
  tools.append(full);model.append(tools,stage);frame.loading='eager';
 }
 const brief=document.createElement('section');brief.className='task-presentation-brief';
 const heading=document.createElement('h2');heading.textContent='简要结论';brief.append(heading);
 const paragraph=document.createElement('p');
 let short=originalReview||originalSummary;
 if(status==='未开始')short=originalSummary||'尚未开始，暂无可交付的模型或训练结果。';
 if(task==='M01')short=[...main.querySelectorAll('#m01-results-overview > .m01-overview-status')].at(-1)?.textContent||'全量TS3 / MGN固定20轮结果已交付，新纯拉伸消融与续训单独记录。';
 if(short.length>230){const sentences=short.match(/[^。！？]+[。！？]?/g)||[short];short='';for(const part of sentences){if(short.length+part.length>230)break;short+=part;}if(!short)short=sentences[0];}
 paragraph.textContent=short||'当前结果、范围和审阅依据保留在下方详情中。';brief.append(paragraph);
 const overview=task==='M01'?main.querySelector(':scope > #m01-results-overview'):null;
 let gallery=overview;
 if(!gallery){
  gallery=document.createElement('section');gallery.className='task-visual-summary';gallery.setAttribute('aria-label','汇总图片');
  const seen=new Set();const preferred=originalImages.filter(img=>/overview|learning|curves|interface|sensitivity|coverage|response|parameters/.test(img.getAttribute('src')));
  const previews=preferred.length?preferred:originalImages;
  for(const img of previews){const src=img.getAttribute('src');if(seen.has(src))continue;seen.add(src);if(seen.size>3)break;
   const figure=document.createElement('figure');figure.className='task-visual-card';const link=document.createElement('a');link.href=src;link.target='_blank';link.rel='noopener';
   const preview=img.cloneNode(false);preview.removeAttribute('id');preview.removeAttribute('style');preview.loading='lazy';link.append(preview);figure.append(link);
   const caption=document.createElement('figcaption');caption.textContent=img.alt||'实际结果汇总图 · 点击查看原图';figure.append(caption);gallery.append(figure);
  }
 }
 const parameterTasks=new Set(['M01','R06','R07','R08']);
 if(parameterTasks.has(task)){
  const figure=document.createElement('figure');figure.className=task==='M01'?'m01-result-card task-parameter-card':'task-visual-card task-parameter-card';
  const link=document.createElement('a');link.href='assets/task-presentation-v1/'+task+'-parameters.svg';link.target='_blank';link.rel='noopener';
  const img=document.createElement('img');img.src=link.href;img.alt=task+' · 数据接口或训练参数图';img.loading='lazy';link.append(img);figure.append(link);
  const caption=document.createElement('figcaption');caption.textContent='配置与训练预算图；原始参数、不同实验口径及来源在详情中。';figure.append(caption);gallery.append(figure);
  const source=document.createElement('p'),sourceLink=document.createElement('a');sourceLink.href='assets/task-presentation-v1/parameters-provenance.json';sourceLink.textContent='参数图的配置来源与校验记录';source.append(sourceLink);content.append(source);
 }
 const keep=new Set([title,stage,overview,detail,...discussions]);
 originalChildren.forEach(node=>{if(!keep.has(node))content.append(node);});
 looseSections.forEach(node=>content.append(node));
 // A secondary details block (e.g. predictions) remains a nested disclosure with all IDs intact.
 main.append(title);if(model)main.append(model);main.append(brief);
 if(gallery.children.length)main.append(gallery);
 else if(!model){const empty=document.createElement('p');empty.className='task-presentation-empty';empty.textContent='暂无已交付模型或汇总图；计划和验收要求见详情。';main.append(empty);}
 main.append(detail,...discussions);main.dataset.presentation='model-first';document.body.classList.add('task-presentation-ready');
 function revealHash(){if(!location.hash)return;const id=decodeURIComponent(location.hash.slice(1));let node=document.getElementById(id);while(node){if(node.tagName==='DETAILS')node.open=true;node=node.parentElement;}}
 revealHash();addEventListener('hashchange',revealHash);
})();
