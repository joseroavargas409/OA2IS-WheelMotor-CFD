# Procedencia de los componentes copiados (third_party/)

Los repositorios originales se usaron **solo como fuentes de referencia de solo lectura**.
Se clonaron en un directorio aparte, se registró su estado (HEAD, `git status`, SHA-256 de
todos los archivos) antes de empezar, y al final se verificó que siguen idénticos
(ver `docs/VALIDACION.md`, sección «Integridad de los repositorios originales»).
Este repositorio no tiene submódulos ni rutas que apunten a ellos: todo lo necesario está copiado aquí.

| Copia en este repositorio | Origen (repo @ commit) | Ruta original | Licencia | Uso |
|---|---|---|---|---|
| `third_party/wheel-motor-inspector/dist/model.js` | `joseroavargas409/wheel-motor-inspector` @ `64bd32e69634f7be45afad43be4bbdbbb5725212` | `dist/model.js` | Sin archivo LICENSE en el repositorio; reutilizado por su autor (José Luis Roa, mismo proyecto OA²IS) | Generador procedural del modelo 3D (cotas, piezas, ejes de giro) |
| `third_party/wheel-motor-inspector/dist/assembly-v2.js` | ídem | `dist/assembly-v2.js` | ídem | Cotas del ensamble rev. 02 (radios, posiciones axiales, layout) |
| `third_party/wheel-motor-inspector/dist/assembly-v3.js` | ídem | `dist/assembly-v3.js` | ídem | Ensamble rev. 03 (piezas, `motions`, relaciones de velocidad) |
| `third_party/wheel-motor-inspector/dist/vendor/three.module.js`, `three.core.js`, `BufferGeometryUtils.js` | ídem | `dist/vendor/` | MIT (three.js authors), ver `LICENSE-three.txt` | Necesario para ejecutar el generador en Node y exportar STL |
| `third_party/wheel-motor-inspector/dist/vendor/LICENSE-three.txt` | ídem | `dist/vendor/LICENSE-three.txt` | MIT | Aviso de licencia de three.js |
| `third_party/komatsu-oil-sensors-v4-1/constants.js` | `joseroavargas409/komatsu-oil-sensors-v4-1` @ `8455afb44dc5156d7fd4e92eb62714397023547e` | `dashboard/constants.js` | MIT (© 2026 OA2IS simulation contributors) | Propiedades del aceite: ρ = 900 kg/m³, ν = 2.64238e-4 m²/s (264 cSt a 60 °C), σ = 0.030 N/m |
| `third_party/komatsu-oil-sensors-v4-1/sensor-model.js` | ídem | `dashboard/sensor-model.js` | MIT | Referencia para la definición de sensores virtuales (alpha, presencia, estado) |
| `third_party/komatsu-oil-sensors-v4-1/LICENSE` | ídem | `LICENSE` | MIT | Texto de la licencia |

Del repositorio `komatsu-oil-sensors-v4` (V4, commit `675d92473857ad0b6bf6dfae5456993fa4cea17e`)
no se copió ningún archivo: solo se inspeccionó (su modelo de fluido es una aproximación por
partículas/campos, no CFD VOF, y V4.1 lo reemplaza).

**No se copiaron** las fotografías, figuras de patente, fotogramas ni documentos de terceros de
`wheel-motor-inspector` (su `THIRD_PARTY_NOTICES.md` aclara que no se concede licencia sobre ellos).

Huellas SHA-256 de las copias (deben coincidir con los archivos originales en los commits indicados):

```
56fb662118c6ee3f1a09de2fec2b1bad2b3e15e8090834362f6baa8ea4f78acb  third_party/komatsu-oil-sensors-v4-1/sensor-model.js
66ab2b90ac214744e77fce2826208c09e700f0aa69617c5fb86573b9c015b7b3  third_party/komatsu-oil-sensors-v4-1/constants.js
2aab283704b028d4505f3709b219d8ef4e6bae19139aa67872f8ae179d0bc3f0  third_party/komatsu-oil-sensors-v4-1/LICENSE
2ed5599e6dc0c140860984467331d2b7ed5042bbd8ea9b7e0d10e6ef7d84c964  third_party/wheel-motor-inspector/dist/model.js
9dc77c856816c9b66272ee2d41eb9532332f02a74a6b62de28175423a83bd951  third_party/wheel-motor-inspector/dist/assembly-v3.js
67d658fac0931a017f8cd198dff7bf0cd905897a5c78a7640e315dfa8ab73902  third_party/wheel-motor-inspector/dist/assembly-v2.js
c8211c69345d2e9949dc7a8ac969380497aa0600a5a8ac6a459c8cd02dd9cb8a  third_party/wheel-motor-inspector/dist/vendor/three.module.js
eb077d2417f61d3e6d9264c317cabc4ea35769ed6b0ab533067292a550784c20  third_party/wheel-motor-inspector/dist/vendor/three.core.js
fda7e946b8e0b5ab39b779206589e7a1079a22eb24efb89d7223e03fdfb1f751  third_party/wheel-motor-inspector/dist/vendor/BufferGeometryUtils.js
bfe119ea4fd413f5f7ca3fcd63adb0c4a073ed39daa2fe7d3e6b769e21272601  third_party/wheel-motor-inspector/dist/vendor/LICENSE-three.txt
```

Las copias no se modificaron. Todo cambio (reparaciones geométricas, regeneración de dientes)
se hace en código nuevo (`oa2is_cfd/model.py`, `oa2is_cfd/gears.py`) y está documentado en
`docs/GEOMETRIA.md`.

## Traducciones (código derivado)

* `oa2is_cfd/oil_sensor.py` es una traducción a Python del modelo capacitivo de
  `third_party/komatsu-oil-sensors-v4-1/sensor-model.js` (MIT, © 2026 OA2IS simulation contributors):
  capacitancia coaxial, película residual, filtro de 0,08 s, ruido xorshift32 con las mismas semillas,
  estados seco/mixto/mojado/no disponible. Verificada contra el original (mismos valores).
* `scenarios/_sensors_v41.yaml` transforma las 12 posiciones de `PROBES` de
  `third_party/komatsu-oil-sensors-v4-1/constants.js` al marco del modelo CFD (`docs/SENSORES.md`).
