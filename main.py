import fdtd.fdtd as fdtd
from fdtd.tools import fft_plot

# grid 100x50===========
sim = fdtd.RCSimulation(
  grid_x=100,
  grid_y=50,
  spatial_diff=1e-3,      # 1mm
  e_r=1,                  # rel. permittivity
  mu_r=1,                 # rel. permeability
  c=3e8,                  # speed of light
  mode=fdtd.RessMode.TM 
)

sim.set_source([fdtd.Sinusoidal(x=50, y=25, freq=10e9)])  # 10GHz
sim.set_probe(75, 25)   # captures field values


# ------ 100x50 - 10GHz - src: 50x25 --------------
# Plane
sim.set_graphics(
  n_steps=700,
  mode=fdtd.GraphicMode.Plane,
  frame_interval=10,
  mpl_cmap="viridis",
  interpol="lanczos",
  title="Resonant Cavity - 100mm x 50mm - 10GHz Sinusoidal Source"
)

sim.render(
  frame_rate=30, 
  dpi=300, 
  output_file="res_cav_100x50_10GHz_50x25_plane.mp4"
)

# Mesh
sim.clear()
sim.set_graphics(
  n_steps=700,
  mode=fdtd.GraphicMode.Mesh,
  frame_interval=30,
  mpl_cmap="viridis",
  title="Resonant Cavity - 100mm x 50mm - 10GHz Sinusoidal Source"
)

sim.render(
  frame_rate=30, 
  dpi=300, 
  output_file="res_cav_100x50_10GHz_50x25_mesh.mp4"
)

fft_plot(sim.get_probe_values(), sim.get_temp_diff(), output_file="fft_res_cav_100x50_10GHz_50x25.png")



# grid 100x100==========
sim = fdtd.RCSimulation(
  grid_x=100,
  grid_y=100,
  spatial_diff=1e-3,      # 1mm
  e_r=1,                  # rel. permittivity
  mu_r=1,                 # rel. permeability
  c=3e8,                  # speed of light
  mode=fdtd.RessMode.TM 
)

sim.set_source([fdtd.Sinusoidal(x=50, y=25, freq=10e9)])  # 10GHz
sim.set_probe(75, 25)   # captures field values


# ------ 100x100 - 10GHz - src: 50x25 --------------
# Plane
sim.set_graphics(
  n_steps=700,
  mode=fdtd.GraphicMode.Plane,
  frame_interval=10,
  mpl_cmap="viridis",
  interpol="lanczos",
  title="Resonant Cavity - 100mm x 100mm - 10GHz Sinusoidal Source"
)

sim.render(
  frame_rate=30, 
  dpi=300, 
  output_file="res_cav_100x100_10GHz_50x25_plane.mp4"
)

# Mesh
sim.clear()
sim.set_graphics(
  n_steps=700,
  mode=fdtd.GraphicMode.Mesh,
  frame_interval=30,
  mpl_cmap="viridis",
  title="Resonant Cavity - 100mm x 100mm - 10GHz Sinusoidal Source"
)

sim.render(
  frame_rate=30, 
  dpi=300, 
  output_file="res_cav_100x100_10GHz_50x25_mesh.mp4"
)

fft_plot(sim.get_probe_values(), sim.get_temp_diff(), output_file="fft_res_cav_100x100_10GHz_50x25.png")


# ------ 100x100 - Gaussian - src: 50x50 --------------
sim.clear()
sim.set_source([fdtd.Gaussian(x=50, y=50, t_center=30, spread=6, amp=2)])

# Plane
sim.set_graphics(
  n_steps=700,
  mode=fdtd.GraphicMode.Plane,
  frame_interval=10,
  mpl_cmap="viridis",
  interpol="lanczos",
  title="Resonant Cavity - 100mm x 100mm - Gaussian Source"
)

sim.render(
  frame_rate=30, 
  dpi=300, 
  output_file="res_cav_100x100_gaussian_50x50_plane.mp4"
)

# Mesh
sim.clear()
sim.set_graphics(
  n_steps=700,
  mode=fdtd.GraphicMode.Mesh,
  frame_interval=30,
  mpl_cmap="viridis",
  title="Resonant Cavity - 100mm x 100mm - Gaussian Source"
)

sim.render(
  frame_rate=30, 
  dpi=300, 
  output_file="res_cav_100x100_gaussian_50x50_mesh.mp4"
)

fft_plot(sim.get_probe_values(), sim.get_temp_diff(), output_file="fft_res_cav_100x100_gaussian_50x50.png")
