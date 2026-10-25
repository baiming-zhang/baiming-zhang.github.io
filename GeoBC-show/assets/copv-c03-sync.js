'use strict';
(() => {
  let pending=null;
  const snapshot=()=>({az,el,roll,zoom,pan:pan.slice(),field,deformation,state:Number(document.getElementById('state').value),cut,wire,bc:document.getElementById('bc-toggle').getAttribute('aria-pressed')==='true',loads:document.getElementById('load-toggle').getAttribute('aria-pressed')==='true'});
  const send=()=>{if(parent!==window)parent.postMessage({kind:'c03-sync',state:snapshot()},location.origin)};
  function apply(s){
    pending=s;if(!manifest)return;az=s.az;el=s.el;roll=s.roll;zoom=s.zoom;pan=s.pan.slice();field=s.field;deformation=s.deformation;
    document.getElementById('field').value=String(field);document.getElementById('deformation').value=String(deformation);document.getElementById('state').value=String(s.state);labels();
    for(const [key,id]of [['bc','bc-toggle'],['loads','load-toggle']]){const b=document.getElementById(id);if(!b.disabled&&(b.getAttribute('aria-pressed')==='true')!==s[key])b.onclick()}
    if(cut!==s.cut)document.getElementById('cut').onclick();if(wire!==s.wire)document.getElementById('mesh-toggle').onclick();draw();
  }
  canvas.addEventListener('pointermove',e=>{if(e.buttons)send()});canvas.addEventListener('wheel',send);canvas.addEventListener('dblclick',send);
  document.querySelector('header').addEventListener('change',send);document.querySelector('header').addEventListener('input',send);document.querySelector('header').addEventListener('click',e=>{if(e.target.closest('button'))send()});
  window.addEventListener('message',e=>{if(e.origin===location.origin&&e.source===parent&&e.data?.kind==='c03-apply')apply(e.data.state)});
  window.addEventListener('c03-ready',()=>{if(pending)apply(pending);parent.postMessage({kind:'c03-ready'},location.origin)});
  window.addEventListener('c03-bc-ready',()=>{if(pending)apply(pending)});
})();
