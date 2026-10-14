import * as THREE from 'three';
import {createModel as createBase,MODEL_RATIO,layout} from './assembly-v2.js';
import {ringGeometry,cylinder} from './model.js';
export {MODEL_RATIO,layout};
const PI=Math.PI,TAU=PI*2;
const yellow='#d9ae37',steel='#b4c0c7',dark='#64747e',violet='#b39bc3',insulation='#dccdae';
// Reference reconstruction. Counts, profiles and fasteners are illustrative, not OEM CAD.
export function createModel(){
 const model=createBase(),{assembly,parts,motions}=model;
 assembly.name='WheelMotor_GDY106_reference_v3';
 Object.assign(assembly.userData,{revision:3,parking_brake:'dry multiple discs; spring applied, hydraulically released',variant:'GDY106 reference; GDY108 catalog only',maintenance_review:'117 sampled visual frames; no full audio transcription'});
 const mat=(color)=>new THREE.MeshStandardMaterial({color,metalness:color===insulation?.05:.72,roughness:color===insulation?.75:.36,side:THREE.DoubleSide});
 function add(p,g,color=p.color,x=0,y=0,z=0,parent=p.group,role=''){const m=new THREE.Mesh(g,mat(color));m.position.set(x,y,z);m.name=p.name+' / '+(p.meshes.length+1);m.userData={part:p.id,baseColor:color,role};parent.add(m);p.meshes.push(m);return m;}
 function ring(p,ro,ri,d,z=0,color=p.color,parent=p.group,role=''){return add(p,ringGeometry(ro,ri,d),color,0,0,z,parent,role);}
 function at(p,g,r,a,z,color=p.color,parent=p.group,role=''){const m=add(p,g,color,r*Math.cos(a),r*Math.sin(a),z,parent,role);m.rotation.z=a;return m;}
 function remove(p,list){for(const m of list){m.parent.remove(m);p.meshes.splice(p.meshes.indexOf(m),1);}}
 function create(id,name,category,color,evidence,description,detail,source,z,explosion,shell=false,visible=true){const group=new THREE.Group();group.name=name;group.position.z=z;assembly.add(group);const p={id,name,category,color,evidence,description,detail,source,group,meshes:[],base:z,explosion,shell,visible,number:parts.size+1};group.userData={part:id,description,detail,source,evidence};parts.set(id,p);return p;}
 function metadata(p,data){Object.assign(p,data);Object.assign(p.group.userData,{description:p.description,detail:p.detail,source:p.source,evidence:p.evidence});}
 function vents(ro,ri,inner,outer,n,width,d){const s=new THREE.Shape();s.absarc(0,0,ro,0,TAU,false);const h=new THREE.Path();h.absarc(0,0,ri,0,TAU,true);s.holes.push(h);for(let i=0;i<n;i++){const a=i*TAU/n,hole=new THREE.Path();hole.moveTo(inner*Math.cos(a-width),inner*Math.sin(a-width));hole.absarc(0,0,inner,a-width,a+width,false);hole.lineTo(outer*Math.cos(a+width),outer*Math.sin(a+width));hole.absarc(0,0,outer,a+width,a-width,true);hole.closePath();s.holes.push(hole);}const g=new THREE.ExtrudeGeometry(s,{depth:d,curveSegments:64,bevelEnabled:true,bevelSize:.002,bevelThickness:.001,bevelSegments:1});g.translate(0,0,-d/2);return g;}
 function shellProfile(points){const g=new THREE.LatheGeometry(points.map(p=>new THREE.Vector2(...p)),96);g.rotateX(PI/2);return g;}
 // Smooth rolled casing and an open, flared connection to the large mounting flange.
 const casing=parts.get('housing');remove(casing,[casing.meshes[0],casing.meshes[1],...casing.meshes.slice(4)]);
 ring(casing,.523,.450,.036,-.355,yellow);
 add(casing,shellProfile([[.481,-.321],[.485,-.33],[.508,-.317],[.521,-.27],[.506,-.21],[.484,-.09],[.484,.30],[.482,.326],[.458,.326],[.457,.30],[.457,-.20],[.480,-.285],[.481,-.321]]));
 add(casing,vents(.486,.12,.244,.414,8,.255,.028),yellow,0,0,.397);
 metadata(casing,{detail:'Cuerpo hueco de perfil suave, testero de ocho ventanas y juntas. Proporciones aproximadas de las fotografías.',source:'i2.png, imagen 930E-4 y ensayos GDY106; perfil reconstruido.',evidence:'Exterior observado; perfil aproximado'});
 const flange=parts.get('flange');remove(flange,flange.meshes.slice(1));
 add(flange,vents(.632,.286,.47,.590,12,.18,.052),yellow,0,0,.02);
 ring(flange,.540,.286,.045,.06);
 for(let i=0;i<12;i++){const a=i*TAU/12+PI/12;const rib=new THREE.CylinderGeometry(.521,.632,.22,96,1,true,a-.024,.048);rib.rotateX(PI/2);add(flange,rib,yellow,0,0,.167);}
 ring(flange,.523,.493,.025,.272);
 metadata(flange,{description:'Brida fija y campana abierta que unen el cuerpo del motor con el soporte. Las ventanas dejan inspeccionar el extremo del conjunto.',detail:'Doce ventanas de sector, costillas de campana y anillo exterior mecanizado con sesenta taladros. Cantidades aproximadas.',source:'Imagen 930E-4, i2.png y CIMETP 00:49; no se dispone de cotas OEM.'});
 // Two bolted joints visible in the supplied exterior photographs.
 const joints=create('jointBands','Bridas y uniones del cubo','Estructura',yellow,'Detalle visual aproximado','Uniones circunferenciales del conjunto exterior. Se modelan como parte del cubo de referencia.','Dos coronas de pernos, arandelas, juntas y banda metálica; 48 pernos por círculo, cantidad ilustrativa.','Imagen 930E-4, i2.png y CIMETP 00:49; pertenencia exacta del montaje no autenticada.',layout.hub,-.75,true);
 const bolt=cylinder(.017,.017,.023,6),washer=ringGeometry(.022,.010,.005);
 for(const worldZ of [.04,-.96]){const z=worldZ-layout.hub;ring(joints,.638,.566,.042,z);ring(joints,.639,.588,.005,z+.024,steel);ring(joints,.592,.566,.034,z-.04,dark);for(let i=0;i<48;i++){const a=i*TAU/48;at(joints,bolt,.611,a,z+.041,steel);at(joints,washer,.611,a,z+.027,steel);}}
 ring(joints,.602,.566,.132,-.88-layout.hub,steel);motions.push({object:joints.group,ratio:1/MODEL_RATIO});
 // Optional windowed OUTER structure seen in CIMETP. The oil enclosure remains sealed.
 const frame=create('outerFrame','Tambor exterior ventilado','Estructura',yellow,'Referencia CIMETP; variante sin identificar','Estructura exterior con ventanas longitudinales vista en el video CIMETP. La envolvente interior del reductor permanece cerrada.','Doce paneles curvos con aberturas reales. Configuración visual opcional; no se atribuye a GDY108 ni a una versión específica.','CIMETP 00:49 y 01:09; no hay placa de serie legible.',-.46,-.75,true,false);
 const R=.600,width=TAU*R/12-.006,L=.91;
 function rounded(path,x,y,w,h,r){path.moveTo(x+r,y);path.lineTo(x+w-r,y);path.quadraticCurveTo(x+w,y,x+w,y+r);path.lineTo(x+w,y+h-r);path.quadraticCurveTo(x+w,y+h,x+w-r,y+h);path.lineTo(x+r,y+h);path.quadraticCurveTo(x,y+h,x,y+h-r);path.lineTo(x,y+r);path.quadraticCurveTo(x,y,x+r,y);path.closePath();}
 const shape=new THREE.Shape();rounded(shape,-width/2,-L/2,width,L,.014);const hole=new THREE.Path();rounded(hole,-.091,-.323,.182,.646,.028);shape.holes.push(hole);
 const panel=new THREE.ExtrudeGeometry(shape,{depth:.018,bevelEnabled:false,curveSegments:12});const pa=panel.attributes.position;for(let i=0;i<pa.count;i++){const u=pa.getX(i),v=pa.getY(i),d=pa.getZ(i),a=u/R;pa.setXYZ(i,(R+d)*Math.cos(a),(R+d)*Math.sin(a),v);}panel.computeVertexNormals();
 for(let i=0;i<12;i++){const mesh=add(frame,panel);mesh.rotation.z=i*TAU/12;}
 for(const z of [-L/2,L/2])ring(frame,.618,.580,.026,z);motions.push({object:frame.group,ratio:1/MODEL_RATIO});
 // Service plugs and lifting points are visual references, not a servicing specification.
 const plugs=create('servicePlugs','Tapones de inspección','Auxiliares',steel,'Detalle ilustrativo','Tapones exteriores para representar accesos de servicio del reductor. Su ubicación exacta exige el manual de la unidad.','Tres alojamientos, tapones hexagonales y arandelas. No se simula nivel ni circulación de aceite.','Fotografías exteriores y referencias de mantenimiento; ubicaciones aproximadas.',layout.hub,-.75);
 for(const [a,z] of [[.7,.42],[3.8,-.28],[5.3,-.80]]){for(const [ro,d,r] of [[.029,.025,.576],[.023,.012,.595]]){const m=at(plugs,cylinder(ro,ro,d,ro===.023?6:24),r,a,z);m.rotation.y=PI/2;}}motions.push({object:plugs.group,ratio:1/MODEL_RATIO});
 const lugs=create('liftingLugs','Orejas de carcasa','Estructura',yellow,'Exterior observado; forma aproximada','Orejas de manipulación visibles en el cuerpo del motor. No se asignan cargas ni procedimientos de izaje.','Dos placas con taladro y nervio de unión.','i2.png e imagen 930E-4 aportadas.',.85,.80);
 for(const a of [-PI/2,PI/2]){const m=at(lugs,ringGeometry(.074,.025,.024),.513,a,-.30);m.rotation.y=PI/2;const n=at(lugs,new THREE.BoxGeometry(.076,.032,.074),.475,a,-.30);n.userData.role='lug_root';}
 // Electrical details derived from the maintenance seminar (embarrilado, 01:30:00).
 const liners=create('slotLiners','Aislación de ranuras','Motor eléctrico',insulation,'Detalle de proceso reconstruido','La aislación separa los conductores del paquete de hierro. Se distinguen los lados rectos y sus salidas hacia las cabezas.','48 revestimientos en U y prolongaciones de salida. Espesor y materiales ilustrativos.','Mantenimiento EVvYz-qH9VE 01:30:00: embarrilado y aislamiento; ranuras del modelo.',.85,1.10);
 const side=new THREE.BoxGeometry(.047,.0014,.558),base=new THREE.BoxGeometry(.0014,.023,.558);
 for(let i=0;i<48;i++){const a=(i+.5)*TAU/48;for(const sign of [-1,1]){const m=at(liners,side,.338,a,0);m.position.x+=-Math.sin(a)*sign*.0115;m.position.y+=Math.cos(a)*sign*.0115;}at(liners,base,.362,a,0);}
 const lacing=create('windingLacing','Amarres de bobinado','Motor eléctrico',insulation,'Proceso visible; patrón aproximado','Bandas de amarre sujetan los grupos de conductores en las cabezas del estator.','24 amarres por extremo; patrón simplificado de cordón continuo. No representa el número OEM de espiras.','Mantenimiento EVvYz-qH9VE 01:30:00: fotografía de embarrilado.',.85,1.30);
 for(const sign of [-1,1])for(let i=0;i<24;i++){const a=(2*i+1)*TAU/48,points=Array.from({length:33},(_,k)=>{const t=k*TAU/32,r=.340+.018*Math.cos(t),z=sign*.320+.017*Math.sin(t);return new THREE.Vector3(r*Math.cos(a),r*Math.sin(a),z);});add(lacing,new THREE.TubeGeometry(new THREE.CatmullRomCurve3(points,true),32,.0018,5,true));}
 const terminals=create('phaseLinks','Puentes de terminales','Motor eléctrico','#cf8854','Detalle de conexión ilustrativo','Tres conexiones de fase unen el estator con la caja de terminales. La ruta es esquemática y no prescribe conexión eléctrica.','Puentes de cobre, soportes aislantes y herrajes de tres fases.','Mantenimiento: cabezas de bobinado; terminales TA/TB/TC de referencias del proyecto.',.85,1.30);
 for(let i=0;i<3;i++){const x=.265+i*.045;add(terminals,new THREE.BoxGeometry(.032,.025,.040),insulation,x,.394,-.005);add(terminals,new THREE.BoxGeometry(.027,.004,.050),'#cf8854',x,.412,-.005);const pts=[new THREE.Vector3(.28+(i-1)*.025,.22,.285),new THREE.Vector3(x,.34,.19),new THREE.Vector3(x,.402,.015)];add(terminals,new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts),32,.005,8,false),'#cf8854');}
 // OEM brochure explicitly specifies spring-applied, hydraulic-release dry multi-disc parking brakes.
 const parking=parts.get('parkingBrake');const oldObjects=new Set();parking.group.traverse(o=>oldObjects.add(o));for(let i=motions.length-1;i>=0;i--)if(oldObjects.has(motions[i].object))motions.splice(i,1);parking.group.clear();parking.meshes=[];
 const rotating=new THREE.Group();rotating.name='Discos secos ligados al rotor';parking.group.add(rotating);motions.push({object:rotating,ratio:1});
 ring(parking,.071,.047,.172,.018,dark,rotating);
 for(let i=0;i<6;i++){const z=-.034+i*.019;ring(parking,.188,.073,.006,z,violet,rotating,'dry_rotating_disc');ring(parking,.198,.072,.004,z+.009,steel,parking.group,'separator');for(let j=0;j<18;j++){const a=j*TAU/18;at(parking,new THREE.BoxGeometry(.011,.010,.006),.071,a,z,dark,rotating);at(parking,new THREE.BoxGeometry(.009,.013,.004),.197,a,z+.009,steel);}}
 metadata(parking,{description:'Freno multidisco seco en el extremo interior del eje del rotor. Los discos giran con el rotor; los separadores reaccionan contra el alojamiento fijo.',detail:'Seis discos móviles y seis separadores, cubo de acoplamiento y pestañas. La cantidad y geometría internas son ilustrativas.',source:'Ficha oficial Komatsu 930E-4SE AESS778-02: multidisco seco, aplicado por resortes y liberado hidráulicamente.',evidence:'Sistema OEM documentado; geometría aproximada'});
 const brakeHousing=create('parkingHousing','Alojamiento del freno seco','Frenos',yellow,'Sistema documentado; forma aproximada','Alojamiento fijo que recibe la reacción del paquete de estacionamiento y contiene su accionamiento.','Campana hueca, placa de reacción, tapa y doce pernos. Ventana central para inspección del acoplamiento.','Ficha oficial 930E-4SE y extremo del motor observado; forma reconstruida.',1.33,2.8,true);
 ring(brakeHousing,.229,.204,.170,.036);ring(brakeHousing,.229,.078,.015,-.064);ring(brakeHousing,.229,.073,.020,.130);
 for(let i=0;i<12;i++)at(brakeHousing,cylinder(.009,.009,.018,6),.214,i*TAU/12,.15,steel);
 // Fixed adapter connects the brake reaction to the motor end-shield, clear of the encoder wheel.
 ring(brakeHousing,.229,.220,.021,-.078);
 const actuator=create('parkingActuator','Pistón y resortes de freno','Frenos',steel,'Principio OEM; disposición ilustrativa','Los resortes aplican el freno. La presión hidráulica desplaza el pistón para liberarlo; aquí se inspecciona el principio sin simular presión.','Pistón anular, dos sellos y doce resortes helicoidales. Cantidades, rigidez y carrera no autenticadas.','Komatsu 930E-4SE: spring applied, hydraulically released dry multiple disc parking brakes.',1.33,2.8);
 ring(actuator,.201,.073,.015,.087,steel,actuator.group,'parking_piston');for(const r of [.077,.197])add(actuator,new THREE.TorusGeometry(r,.0025,6,80),'#283338',0,0,.088);
 const coilPoints=Array.from({length:145},(_,k)=>{const t=k/144;return new THREE.Vector3(.009*Math.cos(t*TAU*7),.009*Math.sin(t*TAU*7),.057*t);});const coil=new THREE.TubeGeometry(new THREE.CatmullRomCurve3(coilPoints),144,.0018,6,false);
 for(let i=0;i<12;i++){const a=i*TAU/12;at(actuator,coil,.169,a,.096,steel,actuator.group,'parking_spring');}
 // Remove previous rotating hardware: old caliper/one-disc coupling is no longer part of this model.
 const fast=parts.get('fasteners');for(const child of [...fast.group.children]){if(child.userData.follows==='flange')for(const mesh of child.children){mesh.position.x*=.616/.568;mesh.position.y*=.616/.568;}if(child.userData.follows==='parkingBrake'){for(let i=motions.length-1;i>=0;i--)if(motions[i].object===child)motions.splice(i,1);remove(fast,fast.meshes.filter(m=>m.parent===child));fast.group.remove(child);}}
 assembly.updateMatrixWorld(true);for(const p of parts.values()){p.group.visible=p.visible;const b=new THREE.Box3().setFromObject(p.group);p.center=b.getCenter(new THREE.Vector3());p.center.z-=p.base;p.dimensions=b.getSize(new THREE.Vector3());p.meshCount=p.meshes.length;}
 model.count=[...parts.values()].reduce((s,p)=>s+p.meshes.length,0);return model;
}
