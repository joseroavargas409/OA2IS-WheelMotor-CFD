# Geometría: del modelo 3D original al dominio CFD

## 1. Fuente

El modelo 3D proviene de `wheel-motor-inspector` @ `64bd32e` (revisión 03 del ensamble,
`dist/assembly-v3.js`, generado proceduralmente con three.js; ver `third_party/PROVENANCE.md`).
No es un CAD de fabricante: es una reconstrucción paramétrica del motor-rueda Komatsu 930E/980E
(GE GDY106/GDY85) hecha a partir de fotos, videos, patentes y documentos (el propio repositorio
lo marca pieza a pieza con un campo `evidence`, p. ej. «Topología GE; dentado ilustrativo»).
Por eso las cotas son aproximadas y algunas piezas se interpenetran; ambas cosas se tratan abajo.

## 2. Extracción (sin modificar el repositorio original)

`geometry/export_parts.mjs` ejecuta en Node 22 el mismo generador (`createModel`) desde la copia en
`third_party/` y escribe:

* `geometry/raw/<pieza>.stl` — 51 STL (41 piezas, 0,2–133 k triángulos), en metros, eje de la rueda = z;
* `geometry/raw/inventory.json` — por pieza: nombre, categoría, evidencia, caja envolvente y
  **cadena de movimiento** (`motion_chain`: razón de velocidad respecto del rotor, fase y pivote).

```
PATH=/opt/node22/bin:$PATH node --import ./geometry/register.mjs geometry/export_parts.mjs geometry/raw
```

Inventario relevante para el sistema de lubricación (razones respecto del rotor = 1):

| Pieza | Movimiento en el modelo | Uso en CFD |
|---|---|---|
| Solar alta (18 d.) + eje del rotor | 1 | cuerpo inmerso móvil (giro) + pared giratoria del eje |
| Satélites alta ×3 (27 d.) | órbita 0,2 + giro propio −0,5333 relativo al porta | 3 cuerpos inmersos móviles (órbita + giro) |
| Porta-satélites alta + eje intermedio | 0,2 | cuerpo inmerso móvil + pared giratoria |
| Corona alta (72 d.) | fija | banda dentada = cuerpo inmerso estático; cuerpo = pared fija |
| Solar baja (23 d.) | 0,2 | cuerpo inmerso móvil |
| Satélites baja ×3 (52 d.) | órbita 0,030667 + giro −0,0749 | 3 cuerpos inmersos móviles |
| Porta-satélites baja = cubo | 0,030667 | cuerpo inmerso + anillo exterior como pared giratoria |
| Corona baja (127 d.) | fija | banda dentada inmersa estática; cuerpo = pared fija |
| Spindle (soporte), asientos | fijos | paredes fijas ajustadas a la malla |
| Cubo de rueda (tapa y tubo) | 0,030667 | pared giratoria (`rotatingWallVelocity`) |

Razón total del modelo: 1/(0,2 · 0,15335) = **32,6087** (igual a la de la reducción real del 930E, 32,62:1).

## 3. Cinemática verificada

`oa2is_cfd/kinematics.py` reproduce la cinemática del original (centros de satélites, fases de
engrane `planet_phase`, giros relativos). `tests/test_kinematics_openfoam.py` la verifica dentro de
OpenFOAM: mueve una malla con `multiSolidBodyMotionSolver` (`multiMotion` = órbita + giro) y
compara los centroides de cada zona con la pose analítica. Resultado: **error máximo 1,5 × 10⁻⁸ m**
tras −2,356 rad del rotor (PASSED).

## 4. Del modelo visual a un dominio de fluido válido

El dominio del aceite es la cavidad de los engranajes: radialmente hasta la cara interior del cubo
(r = 0,536 m; 0,516 m en la brida de la tapa) y axialmente desde la cara interior de la tapa del cubo
(z = −1,7345 m) hasta la cara exterior del asiento del rodamiento del cubo (z = −0,8925 m).
`oa2is_cfd/model.py` describe cada sólido como un prisma 2D extruido axialmente (todas las piezas
del original son extrusiones o revoluciones), con la cota de origen anotada en cada uno.

Reparaciones necesarias (lista `REPAIRS` en `oa2is_cfd/model.py`):

1. **Dientes regenerados como evolventes exactas con juego (backlash)**. En el original los perfiles
   de la corona están espejados y los dientes se interpenetran (juego cero): un fluido no puede existir
   entre dientes que se solapan. `oa2is_cfd/gears.py` genera solar/satélite/corona con el mismo número
   de dientes, radio primitivo y fase de engrane, con juego 4 mm (alta) y 3 mm (baja); se verificó
   solape cero en un ciclo completo de engrane (holgura mínima 1,88 / 1,41 mm).
2. Soporte anular del spindle en el lado exterior de la etapa baja recortado a [−1,559, −1,548]:
   en el original se superpone con el anillo giratorio del porta-satélites bajo.
3. Arandelas de empuje en z ≈ −1,22 omitidas: el original las deja fijas mientras intersectan
   la placa giratoria del porta alto y la brida del eje intermedio.
4. Interior del manguito intermedio (r < 0,049) y perforaciones del eje/spindle tratados como
   volúmenes cerrados (sellados).
5. Radios a menos de 1–3 mm entre sí ajustados para poder mallar (p. ej. tubo del spindle
   0,532 → 0,536, que cierra una rendija de 4 mm entre tubo y cubo); detalle en `REPAIRS`.
6. Cavidad cerrada en la cara interior de la tapa y en la cara exterior del asiento del rodamiento;
   lo que queda más allá pertenece a la tapa de empuje y a la cavidad del rodamiento.

La geometría procesada (pose de reposo) se exporta con

```
python3 tools/export_cfd_geometry.py          # geometry/cfd/*.stl + geometry_summary.json
xvfb-run -a python3 tools/render_geometry.py docs/img/geometria_cfd.png
```

![geometría CFD](img/geometria_cfd.png)

## 5. Volumen y nivel de 95 L

* Volumen geométrico de fluido de la cavidad (fuera de todos los sólidos, en reposo): **493,1 L**.
  95 L corresponden a una superficie libre plana en y = −0,272 m (19 % del volumen), por debajo del
  eje: es geométricamente consistente con un cárter húmedo de reductor planetario.
* En la malla de frontera inmersa, las bolsas de fluido que quedan **completamente encerradas**
  entre piezas fijas (65 L, sobre todo la bolsa interior junto al asiento del rodamiento) se retiran
  porque no tienen comunicación con el resto; los 95 L se distribuyen en los 426–428 L comunicados,
  lo que sube el nivel inicial a y = −0,234 m. El nivel se calcula por bisección **sobre la malla
  discreta**, de modo que el aceite inicial integrado es exactamente 95,000 L.
* El escenario de bajo inventario usa 55 L.

## 6. Limitaciones geométricas

* Cotas aproximadas (modelo reconstruido, no CAD OEM). Ancho de cara de engranajes 102 mm (original).
* Rodamientos de satélites, sellos y arandelas se representan como núcleos sólidos o se omiten.
* No hay pasajes internos de lubricación (no existen en el modelo original).
