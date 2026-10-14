// The only experimental factors are oil liters and number of virtual probes.
export const FIXED=Object.freeze({version:'4.1',maturity:'Simulada – no validada',model:'930E-4SE',duration_s:20,frame_interval_s:.02,settle_s:3,ramp_s:2,
  motor_rpm:300,modeled_ratio:750/23,published_ratio:32.62,temperature_c:60,
  gravity:[0,0,-9.81],pitch_deg:0,roll_deg:0,acceleration_m_s2:0,
  density_kg_m3:900,kinematic_viscosity_m2_s:.000264238,particle_radius_m:.016,surface_tension_n_m:.030,wall_wetting_gain_m_s2:.3,
  pressure_solver:'DFSPH',sensor_tau_s:.08,sensor_noise_pf:.01,sensor_seed:930,
  air_model:'SPH viewer infers gas in fluid space; a separately labelled OpenFOAM oil-air reference is provided',
  mechanism:'two serial planetary stages with stationary rings; prescribed constant motor speed',
  geometry_origin:'illustrative dimensions and tooth counts; no OEM CAD or calibrated contact model',
  stage_teeth:[[18,27,72],[23,52,127]],numerical_projection_limit_m:.008});
export const PROBES=Object.freeze([
  {id:'S1',name:'Retorno de alta',stage:'Alta',position:[-.24,.33,-.24]},
  {id:'S2',name:'Lateral de baja',stage:'Baja',position:[.335,-.33,-.16]},
  {id:'S3',name:'Fondo del cárter',stage:'Común',position:[0,0,-.405]},
  {id:'S4',name:'Espacio superior',stage:'Común',position:[0,0,.22]},
  {id:'S5',name:'Lateral de alta',stage:'Alta',position:[.35,.33,-.14]},
  {id:'S6',name:'Retorno de baja',stage:'Baja',position:[-.335,-.33,-.16]},
  {id:'S7',name:'Canal central inferior',stage:'Común',position:[-.12,0,-.28]},
  {id:'S8',name:'Canal central superior',stage:'Común',position:[.13,0,.10]},
  {id:'S9',name:'Zona de entrada',stage:'Alta',position:[0,.38,-.25]},
  {id:'S10',name:'Zona de salida',stage:'Baja',position:[0,-.38,-.25]},
  {id:'S11',name:'Salpicadura superior derecha',stage:'Común',position:[.21,0,.20]},
  {id:'S12',name:'Salpicadura superior izquierda',stage:'Común',position:[-.21,0,.20]}
]);
export const SENSOR_COUNTS=[2,4,6,8,12];
export const VOLUMES=[55,95,125];
export const COLORS=['#5ce0c5','#ff9a69','#8cbbff','#e4a0f8','#e4d673','#86dd92','#ff89b2','#78d8f4','#c4a181','#afa5ff','#d0e19a','#80b5b5'];
