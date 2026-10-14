import {FIXED,PROBES} from './constants.js';
export function kernel(r,h){const q=r/h;if(q>=1)return 0;return 8/(Math.PI*h**3)*(q<=.5?1+6*(q**3-q**2):2*(1-q)**3);}
const clamp=x=>Math.min(1,Math.max(0,x));
export class Geometry {
 constructor(bodies=[]){this.bodies=bodies;}
 nearSurface(q,h,t){
  if(.5-Math.hypot(q[0],q[2])<h||.45-Math.abs(q[1])<h)return true;
  for(const b of this.bodies){if(Math.abs(q[1]-b.center_y)>b.length/2+h)continue;const a=b.phase-b.orbit_factor*FIXED.motor_rpm*2*Math.PI/60*t,r=Math.hypot(q[0]-b.orbit*Math.cos(a),q[2]-b.orbit*Math.sin(a));if(r>b.inner-h&&r<b.outer+h)return true;}return false;
 }
 solid(q,t){
  if(Math.hypot(q[0],q[2])>=.5||Math.abs(q[1])>=.45)return true;
  for(const b of this.bodies){if(Math.abs(q[1]-b.center_y)>=b.length/2)continue;
   const orbitAngle=b.phase-b.orbit_factor*FIXED.motor_rpm*2*Math.PI/60*t;
   const cx=b.orbit*Math.cos(orbitAngle),cz=b.orbit*Math.sin(orbitAngle),dx=q[0]-cx,dz=q[2]-cz,r=Math.hypot(dx,dz);
   if(r<=b.inner||r>=b.outer)continue;
   if(!b.teeth)return true;
   const spin=b.speed_factor*FIXED.motor_rpm*2*Math.PI/60*t;
   const x=Math.cos(spin)*dx-Math.sin(spin)*dz,z=Math.sin(spin)*dx+Math.cos(spin)*dz;
   const angle=((Math.atan2(z,x)%(2*Math.PI))+2*Math.PI)%(2*Math.PI),n=b.teeth*4,scaled=angle*n/(2*Math.PI),j=Math.floor(scaled),f=scaled-j,da=2*Math.PI/n;
   const vertex=i=>b.inward?[b.inner+(i%4>=2?.014:0),b.outer]:[b.inner,b.outer-(i%4>=2?.014:0)];
   const [a,c]=vertex(j),[d,e]=vertex((j+1)%n),beta=f*da;
   const ray=(r1,r2)=>r1*r2*Math.sin(da)/(r2*Math.sin(da-beta)+r1*Math.sin(beta));
   if(r>ray(a,d)&&r<ray(c,e))return true;
  }return false;
 }
}
export function coaxCapacitance(oilFraction,filmThickness=0){
 const eps0=8.8541878128e-12,a=.004,b=.006,length=.025,epsOil=2.2,epsGas=1.,stray=8e-12;
 const film=Math.max(0,Math.min(.00015,filmThickness));
 const mixture=epsGas+(epsOil-epsGas)*clamp(oilFraction);
 const inverse=Math.log((a+film)/a)/epsOil+Math.log((b-film)/(a+film))/mixture+Math.log(b/(b-film))/epsOil;
 return (stray+2*Math.PI*eps0*length/inverse)*1e12;
}
const OFFSETS=[[0,0,0],[.006,0,0],[-.006,0,0],[0,.006,0],[0,-.006,0],[0,0,.0125],[0,0,-.0125],[0,0,.00625],[0,0,-.00625]];
export class SensorArray {
 constructor(particleVolume,{particleRadius=FIXED.particle_radius_m,geometry=[]}={}){
  this.volume=particleVolume;this.radius=particleRadius;this.h=4*particleRadius;this.geometry=new Geometry(geometry);this.lastTime=null;
  this.filtered=new Float64Array(PROBES.length);this.film=new Float64Array(PROBES.length);this.initialized=false;
  this.seeds=PROBES.map((_,i)=>(FIXED.sensor_seed+(i+1)*2654435761)>>>0);
  const spacing=2*particleRadius;let norm=0;for(let x=-3;x<=3;x++)for(let y=-3;y<=3;y++)for(let z=-3;z<=3;z++)norm+=kernel(Math.hypot(x,y,z)*spacing,this.h)*particleVolume;
  this.norm=norm;this.stencil=[];for(let x=-2;x<=2;x++)for(let y=-2;y<=2;y++)for(let z=-2;z<=2;z++){const dx=x*spacing,dy=y*spacing,dz=z*spacing,w=kernel(Math.hypot(dx,dy,dz),this.h)*particleVolume;if(w>0)this.stencil.push([dx,dy,dz,w]);}
 }
 random(i){let a=this.seeds[i];a^=a<<13;a^=a>>>17;a^=a<<5;this.seeds[i]=a>>>0;return ((a>>>0)+.5)/4294967296;}
 sample(time,positions,mechanismTime=time){
  const h=this.h,grid=new Map(),key=(x,y,z)=>`${x},${y},${z}`;
  for(let j=0;j<positions.length;j+=3){const cell=key(Math.floor(positions[j]/h),Math.floor(positions[j+1]/h),Math.floor(positions[j+2]/h));let entry=grid.get(cell);if(!entry){entry={ids:[],x:Math.floor(positions[j]/h),y:Math.floor(positions[j+1]/h),z:Math.floor(positions[j+2]/h),connected:false};grid.set(cell,entry);}entry.ids.push(j);}
  const pending=[];for(const e of grid.values())if(e.z*h<-.42&&Math.abs(e.x*h)<.3){e.connected=true;pending.push(e);}
  const neighbors=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]];
  for(let i=0;i<pending.length;i++){const e=pending[i];for(const d of neighbors){const f=grid.get(key(e.x+d[0],e.y+d[1],e.z+d[2]));if(f&&!f.connected){f.connected=true;pending.push(f);}}}
  const occupancy=q=>{let total=0,connected=0;if(this.geometry.solid(q,mechanismTime))return {oil:0,bulk:0,solid:true};
   const base=q.map(v=>Math.floor(v/h));for(let x=-1;x<=1;x++)for(let y=-1;y<=1;y++)for(let z=-1;z<=1;z++){const e=grid.get(key(base[0]+x,base[1]+y,base[2]+z));if(!e)continue;for(const j of e.ids){const weight=kernel(Math.hypot(positions[j]-q[0],positions[j+1]-q[1],positions[j+2]-q[2]),h)*this.volume;total+=weight;if(e.connected)connected+=weight;}}
   let norm=this.norm;if(this.geometry.nearSurface(q,h,mechanismTime)){norm=0;for(const [x,y,z,w] of this.stencil)if(!this.geometry.solid([q[0]+x,q[1]+y,q[2]+z],mechanismTime))norm+=w;}
   return {oil:clamp(total/Math.max(norm,1e-8)),bulk:clamp(connected/Math.max(norm,1e-8)),solid:false};
  };
  const dt=this.lastTime===null?0:Math.max(0,time-this.lastTime),alpha=this.lastTime===null?1:1-Math.exp(-dt/FIXED.sensor_tau_s);this.lastTime=time;
  const result=PROBES.map((s,i)=>{let oil=0,bulk=0,valid=0;
   for(const d of OFFSETS){const v=occupancy(s.position.map((q,k)=>q+d[k]));if(!v.solid){oil+=v.oil;bulk+=v.bulk;valid++;}}
   if(!valid)return {id:s.id,oil_fraction:null,air_fraction:null,bulk_fraction:null,detached_fraction:null,film_thickness_um:null,solid_excluded_fraction:1,signal:null,capacitance_proxy_pf:null,capacitance_pf:null,ideal_capacitance_pf:null,level_mm:null,equivalent_oil_column_mm:null,wet:false,state:'unavailable',temperature_c:FIXED.temperature_c};
   oil/=valid;bulk/=valid;
   const targetFilm=.00005*oil,filmTau=targetFilm>this.film[i]?.1:.7;
   if(!this.initialized)this.film[i]=targetFilm;else this.film[i]+=(targetFilm-this.film[i])*(1-Math.exp(-dt/filmTau));
   const target=coaxCapacitance(oil,this.film[i]);if(!this.initialized)this.filtered[i]=target;else this.filtered[i]+=alpha*(target-this.filtered[i]);
   const noise=Math.sqrt(-2*Math.log(this.random(i)))*Math.cos(2*Math.PI*this.random(i))*.01;
   const measured=this.filtered[i]+noise,dry=coaxCapacitance(0),wet=coaxCapacitance(1),signal=clamp((measured-dry)/(wet-dry));
   let bulkLevel=0,equivalentLevel=0;const spacing=2*this.radius;
   for(let z=-.5+this.radius;z<=.5-this.radius;z+=spacing){const value=occupancy([s.position[0],s.position[1],z]);if(!value.solid){bulkLevel+=value.bulk*spacing*1000;equivalentLevel+=value.oil*spacing*1000;}}
   return {id:s.id,oil_fraction:oil,air_fraction:1-oil,air_basis:'inferred_from_liquid_field_in_fluid_space',
    bulk_fraction:bulk,detached_fraction:Math.max(0,oil-bulk),film_thickness_um:this.film[i]*1e6,solid_excluded_fraction:1-valid/OFFSETS.length,
    signal,capacitance_proxy_pf:measured,capacitance_pf:measured,ideal_capacitance_pf:target,
    level_mm:bulkLevel,equivalent_oil_column_mm:equivalentLevel,wet:signal>.6,state:!valid?'unavailable':signal<.2?'dry':signal<.8?'mixed':'wet',temperature_c:FIXED.temperature_c};
  });this.initialized=true;return result;
 }
}
