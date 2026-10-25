'use strict';
// Closed triangle meshes, instanced in model coordinates. No arrow sprites,
// screen-facing quads, SVG arrowheads, or outline pass.
window.FESpatialGlyphs=class {
  constructor(gl){
    this.gl=gl;this.shapes=new Map();
    const shader=(type,source)=>{const s=gl.createShader(type);gl.shaderSource(s,source);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(s));return s;};
    this.program=gl.createProgram();
    gl.attachShader(this.program,shader(gl.VERTEX_SHADER,`#version 300 es
    in vec3 aLocal,aNormal,aAnchor,aU,aDirection,aColor;in vec2 aSize;
    uniform vec3 uCenter,uRight,uUp,uDepth;uniform vec2 uPan;
    uniform float uScale,uAspect,uZScale,uDeform;
    out vec3 normalView,baseColor;
    void main(){vec3 d=normalize(aDirection),helper=abs(d.z)>.8?vec3(0,1,0):vec3(0,0,1);
      vec3 x=normalize(cross(helper,d)),y=cross(d,x);
      vec3 p=aAnchor+aU*uDeform+(x*aLocal.x+y*aLocal.y+d*(aLocal.z+aSize.y))*aSize.x-uCenter;
      vec3 n=x*aNormal.x+y*aNormal.y+d*aNormal.z;
      normalView=vec3(dot(n,uRight),dot(n,uUp),dot(n,uDepth));baseColor=aColor;
      gl_Position=vec4(dot(p,uRight)/(uScale*uAspect)+uPan.x,dot(p,uUp)/uScale+uPan.y,-dot(p,uDepth)/uZScale,1.);
    }`));
    gl.attachShader(this.program,shader(gl.FRAGMENT_SHADER,`#version 300 es
    precision highp float;in vec3 normalView,baseColor;out vec4 color;
    void main(){vec3 n=normalize(normalView),l=normalize(vec3(-.35,.55,1.));float diffuse=max(dot(n,l),0.);
      float shine=pow(max(dot(reflect(-l,n),vec3(0,0,1)),0.),28.);
      color=vec4(baseColor*(.52+.48*diffuse)+vec3(.22)*shine,1.);
    }`));
    gl.linkProgram(this.program);if(!gl.getProgramParameter(this.program,gl.LINK_STATUS))throw Error(gl.getProgramInfoLog(this.program));
    this.instances=gl.createBuffer();
  }
  geometry(kind){
    if(this.shapes.has(kind))return this.shapes.get(kind);
    const g=this.gl,vertices=[],n=16;
    const vertex=(p,normal)=>vertices.push(...p,...normal);
    const triangle=(a,b,c,na,nb=na,nc=na)=>{vertex(a,na);vertex(b,nb);vertex(c,nc);};
    const disk=(z,r,sign)=>{for(let i=0;i<n;i++){const a=i*2*Math.PI/n,b=(i+1)*2*Math.PI/n;triangle([0,0,z],[r*Math.cos(a),r*Math.sin(a),z],[r*Math.cos(b),r*Math.sin(b),z],[0,0,sign]);}};
    const cylinder=(z0,z1,r)=>{for(let i=0;i<n;i++){const a=i*2*Math.PI/n,b=(i+1)*2*Math.PI/n,na=[Math.cos(a),Math.sin(a),0],nb=[Math.cos(b),Math.sin(b),0],p=[r*na[0],r*na[1],z0],q=[r*nb[0],r*nb[1],z0],s=[r*na[0],r*na[1],z1],t=[r*nb[0],r*nb[1],z1];triangle(p,q,s,na,nb,na);triangle(s,q,t,na,nb,nb);}};
    if(kind==='ball'){
      const point=(a,b)=>[Math.sin(a)*Math.cos(b),Math.sin(a)*Math.sin(b),Math.cos(a)];
      for(let j=0;j<8;j++)for(let i=0;i<n;i++){const a=point(j*Math.PI/8,i*2*Math.PI/n),b=point(j*Math.PI/8,(i+1)*2*Math.PI/n),c=point((j+1)*Math.PI/8,i*2*Math.PI/n),d=point((j+1)*Math.PI/8,(i+1)*2*Math.PI/n);triangle(a,b,c,a,b,c);triangle(c,b,d,c,b,d);}
    }else if(kind==='support'){
      cylinder(0,.88,.05);disk(0,.05,-1);cylinder(.88,1,.14);disk(.88,.14,-1);disk(1,.14,1);
    }else{
      cylinder(0,.72,.035);disk(0,.035,-1);disk(.72,.115,-1);
      for(let i=0;i<n;i++){const a=i*2*Math.PI/n,b=(i+1)*2*Math.PI/n,k=.115/.28,normal=t=>{const m=Math.hypot(1,k);return [Math.cos(t)/m,Math.sin(t)/m,k/m];};triangle([.115*Math.cos(a),.115*Math.sin(a),.72],[.115*Math.cos(b),.115*Math.sin(b),.72],[0,0,1],normal(a),normal(b),normal((a+b)/2));}
    }
    const previous=g.getParameter(g.VERTEX_ARRAY_BINDING),vao=g.createVertexArray();g.bindVertexArray(vao);
    const buffer=g.createBuffer();g.bindBuffer(g.ARRAY_BUFFER,buffer);g.bufferData(g.ARRAY_BUFFER,new Float32Array(vertices),g.STATIC_DRAW);
    for(const [name,offset]of [['aLocal',0],['aNormal',12]]){const loc=g.getAttribLocation(this.program,name);g.enableVertexAttribArray(loc);g.vertexAttribPointer(loc,3,g.FLOAT,false,24,offset);}
    g.bindBuffer(g.ARRAY_BUFFER,this.instances);
    for(const [name,size,offset]of [['aAnchor',3,0],['aU',3,12],['aDirection',3,24],['aColor',3,36],['aSize',2,48]]){const loc=g.getAttribLocation(this.program,name);g.enableVertexAttribArray(loc);g.vertexAttribPointer(loc,size,g.FLOAT,false,56,offset);g.vertexAttribDivisor(loc,1);}
    g.bindVertexArray(previous);const result={vao,count:vertices.length/6};this.shapes.set(kind,result);return result;
  }
  draw(frame,items){
    if(!items.length)return;
    const g=this.gl,p=this.program,previous={program:g.getParameter(g.CURRENT_PROGRAM),vao:g.getParameter(g.VERTEX_ARRAY_BINDING),viewport:g.getParameter(g.VIEWPORT),scissor:g.getParameter(g.SCISSOR_BOX),scissorOn:g.isEnabled(g.SCISSOR_TEST),depth:g.isEnabled(g.DEPTH_TEST),blend:g.isEnabled(g.BLEND),cull:g.isEnabled(g.CULL_FACE),depthMask:g.getParameter(g.DEPTH_WRITEMASK),depthFunc:g.getParameter(g.DEPTH_FUNC)};
    const v=frame.viewport;g.viewport(v.x,v.y,v.w,v.h);g.enable(g.SCISSOR_TEST);g.scissor(v.x,v.y,v.w,v.h);
    // X-ray visibility uses a fresh depth layer, so the actual 3D glyphs still
    // occlude one another correctly while internal load locations remain visible.
    g.depthMask(true);g.clear(g.DEPTH_BUFFER_BIT);g.enable(g.DEPTH_TEST);g.depthFunc(g.LEQUAL);g.disable(g.BLEND);g.disable(g.CULL_FACE);g.useProgram(p);
    for(const name of ['center','right','up','depth'])g.uniform3fv(g.getUniformLocation(p,'u'+name[0].toUpperCase()+name.slice(1)),frame[name]);
    g.uniform2fv(g.getUniformLocation(p,'uPan'),frame.pan);
    for(const [name,value]of [['uScale',frame.scale],['uAspect',v.w/v.h],['uZScale',frame.depthScale||frame.radius*4],['uDeform',frame.deform||0]])g.uniform1f(g.getUniformLocation(p,name),value);
    for(const kind of ['arrow','support','ball']){
      const selected=items.filter(x=>(x.shape||'arrow')===kind);if(!selected.length)continue;
      const mesh=this.geometry(kind),data=new Float32Array(selected.length*14);
      selected.forEach((item,i)=>data.set([...item.p,...(item.u||[0,0,0]),...item.d,...item.color,item.length,item.shift||0],i*14));
      g.bindVertexArray(mesh.vao);g.bindBuffer(g.ARRAY_BUFFER,this.instances);g.bufferData(g.ARRAY_BUFFER,data,g.STREAM_DRAW);g.drawArraysInstanced(g.TRIANGLES,0,mesh.count,selected.length);
    }
    g.bindVertexArray(previous.vao);g.useProgram(previous.program);g.viewport(...previous.viewport);g.scissor(...previous.scissor);
    for(const [flag,on]of [[g.SCISSOR_TEST,previous.scissorOn],[g.DEPTH_TEST,previous.depth],[g.BLEND,previous.blend],[g.CULL_FACE,previous.cull]])on?g.enable(flag):g.disable(flag);
    g.depthMask(previous.depthMask);g.depthFunc(previous.depthFunc);
  }
};
