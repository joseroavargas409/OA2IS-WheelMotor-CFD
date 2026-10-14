// Node ESM resolve hook: maps the bare specifier 'three' used by the
// wheel-motor-inspector sources to the vendored Three.js copy, so the copied
// sources stay byte-identical to the originals (see third_party/PROVENANCE.md).
import {pathToFileURL} from 'node:url';
import path from 'node:path';
const here = path.dirname(new URL(import.meta.url).pathname);
const THREE = pathToFileURL(path.resolve(here, '../third_party/wheel-motor-inspector/dist/vendor/three.module.js')).href;
export async function resolve(specifier, context, next) {
  if (specifier === 'three') return {url: THREE, shortCircuit: true};
  return next(specifier, context);
}
