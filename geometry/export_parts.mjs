// Export the wheel-motor-inspector (revision 03) assembly, unmodified, to
// per-part binary STL files in world coordinates (metres, rotation axis = +Z)
// plus a JSON inventory with kinematics (motion ratio of every mesh relative to
// the rotor, and the pivot of orbiting planets).
//
// Usage: node --import ./geometry/register.mjs geometry/export_parts.mjs <outDir>
import fs from 'node:fs';
import path from 'node:path';
import * as THREE from 'three';
import {createModel} from '../third_party/wheel-motor-inspector/dist/assembly-v3.js';
import {stages, layout, MODEL_RATIO} from '../third_party/wheel-motor-inspector/dist/assembly-v2.js';

const out = process.argv[2] || 'geometry/raw';
fs.mkdirSync(out, {recursive: true});
const model = createModel();
model.assembly.updateMatrixWorld(true);

// Map every Object3D that carries a motion to its ratio (rotation about the
// object's own local Z, relative to its parent).
const motionOf = new Map(model.motions.map(m => [m.object, m]));

function chainFor(obj) {
  // List of motion-carrying ancestors from root to obj (inclusive).
  const chain = [];
  for (let o = obj; o; o = o.parent) if (motionOf.has(o)) chain.unshift(o);
  return chain.map(o => {
    const m = motionOf.get(o);
    const p = new THREE.Vector3(); o.getWorldPosition(p);
    return {name: o.name || o.parent?.name || '', ratio: m.ratio, phase: m.phase || 0,
            pivot_world: [p.x, p.y, p.z]};
  });
}

function writeSTL(file, triangles) {
  const buf = Buffer.alloc(84 + 50 * triangles.length);
  buf.write('wheel-motor-inspector rev03 export (metres)', 0, 'ascii');
  buf.writeUInt32LE(triangles.length, 80);
  let o = 84;
  const a = new THREE.Vector3(), b = new THREE.Vector3(), c = new THREE.Vector3();
  for (const t of triangles) {
    a.fromArray(t, 0); b.fromArray(t, 3); c.fromArray(t, 6);
    const n = new THREE.Vector3().subVectors(c, b).cross(new THREE.Vector3().subVectors(a, b)).normalize();
    for (const v of [n, a, b, c]) { buf.writeFloatLE(v.x, o); buf.writeFloatLE(v.y, o + 4); buf.writeFloatLE(v.z, o + 8); o += 12; }
    buf.writeUInt16LE(0, o); o += 2;
  }
  fs.writeFileSync(file, buf);
}

function meshTriangles(mesh) {
  const g = mesh.geometry, pos = g.attributes.position, idx = g.index;
  const M = mesh.matrixWorld, tris = [], v = new THREE.Vector3();
  const n = idx ? idx.count : pos.count;
  for (let i = 0; i < n; i += 3) {
    const t = [];
    for (let k = 0; k < 3; k++) {
      const j = idx ? idx.getX(i + k) : i + k;
      v.fromBufferAttribute(pos, j).applyMatrix4(M); t.push(v.x, v.y, v.z);
    }
    tris.push(t);
  }
  return tris;
}

const inventory = {source: 'wheel-motor-inspector@64bd32e (dist/assembly-v3.js createModel)',
  units: 'm', axis: [0, 0, 1], model_ratio: MODEL_RATIO, layout, stages, parts: []};
for (const p of model.parts.values()) {
  const meshes = [];
  // Group the part's meshes by kinematic chain so orbiting planets get one STL each.
  const groups = new Map();
  for (const mesh of p.meshes) {
    const chain = chainFor(mesh);
    const key = JSON.stringify(chain.map(c => [c.ratio, c.pivot_world.map(x => +x.toFixed(5))]));
    if (!groups.has(key)) groups.set(key, {chain, meshes: []});
    groups.get(key).meshes.push(mesh);
  }
  let gi = 0;
  for (const {chain, meshes: ms} of groups.values()) {
    const tris = ms.flatMap(meshTriangles);
    const file = `${p.id}${groups.size > 1 ? '_' + gi : ''}.stl`;
    const box = new THREE.Box3(); for (const m of ms) box.expandByObject(m);
    writeSTL(path.join(out, file), tris);
    meshes.push({file, mesh_count: ms.length, triangles: tris.length,
      roles: [...new Set(ms.map(m => m.userData.role).filter(Boolean))],
      bbox_min: box.min.toArray(), bbox_max: box.max.toArray(), motion_chain: chain});
    gi++;
  }
  inventory.parts.push({id: p.id, name: p.name, category: p.category, visible_by_default: p.visible,
    evidence: p.evidence, groups: meshes});
}
fs.writeFileSync(path.join(out, 'inventory.json'), JSON.stringify(inventory, null, 1));
console.log(JSON.stringify({parts: inventory.parts.length, files: inventory.parts.reduce((s, p) => s + p.groups.length, 0)}));
