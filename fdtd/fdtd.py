import numpy as np
from enum import Enum
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ------------------ Sources ---------------------------
class Source:
  def __init__(self, x: int, y: int):
    self.x = x
    self.y = y

  def get_state(self, step):
    raise NotImplementedError

class Gaussian(Source):
  def __init__(self, x: int, y: int, spread: int, t_center: int, amp=1):
    super().__init__(x, y)
    self.__spread = spread
    self.__tcenter = t_center
    self.__amp = amp

  def get_state(self, step):
    return self.__amp * np.exp(-0.5 * ((step - self.__tcenter) / self.__spread)**2)
  

class Sinusoidal(Source):
  def __init__(self, x: int, y: int ,freq: float, amp=1):
    super().__init__(x, y)
    self.__freq = freq
    self.__amp = amp

  def get_state(self, step):
    return  self.__amp * np.sin(2 * np.pi * self.__freq * step)


# ------------ Modes --------------------------------------
class RessMode(Enum):
  TM = "TM"
  TE = "TE"

class GraphicMode(Enum):
  Mesh = "mesh"
  Plane = "plane"


# ------------- Simulation ------------------------------
class RCSimulation:
  def __init__(self, grid_x, grid_y, spatial_diff=1e-3, e_r=1, mu_r=1, c=3e8, mode=RessMode.TM):
    self.__epsilon = 8.854187817e-12 * e_r
    self.__mu = 4 * np.pi * 1e-7 * mu_r

    # spatial and temporal parameters
    self.__dx = spatial_diff  
    self.__dy = spatial_diff 
    self.__dt = spatial_diff / (2 * c)  # respecting the CFL condition     
    
    # grid dimensions
    self.__grid_x = grid_x
    self.__grid_y = grid_y 
    self.__n_steps = 0

    # initialize fields
    self.__Ez = np.zeros((grid_x, grid_y))  
    self.__Hx = np.zeros((grid_x, grid_y))  
    self.__Hy = np.zeros((grid_x, grid_y))
    self.__mode = mode

    # for FFT analysis
    self.__ez_history = []  
    self.__probe = (0,0)

    # sources
    self.__sources = []

    # simulation
    self.__anim = None


  # update magnetic fields Hx and Hy based on the current electric field Ez.
  def __update_H(self):
    if self.__mode == RessMode.TM:
      self.__Hx[:, :-1] -= (self.__dt / (self.__mu * self.__dy)) * (self.__Ez[:, 1:] - self.__Ez[:, :-1])
      self.__Hy[:-1, :] += (self.__dt / (self.__mu * self.__dx)) * (self.__Ez[1:, :] - self.__Ez[:-1, :])
    
    elif self.__mode == RessMode.TE:
      print("TE mode not implemented yet.")
      exit(1)


  # update electric field Ez based on the current magnetic fields Hx and Hy.
  def __update_E(self):
    if self.__mode == RessMode.TM:
      self.__Ez[1:-1, 1:-1] += (self.__dt / (self.__epsilon * self.__dx)) * (
          (self.__Hy[1:-1, 1:-1] - self.__Hy[0:-2, 1:-1]) -
          (self.__Hx[1:-1, 1:-1] - self.__Hx[1:-1, 0:-2])
      )
    elif self.__mode == RessMode.TE:
      print("TE mode not implemented yet.")
      exit(1)


  def __apply_source(self, t):
    for src in self.__sources:
      if isinstance(src, Sinusoidal):
        t = t * self.__dt
      self.__Ez[src.x, src.y] += src.get_state(t)


  # apply boundary conditions to the electric field Ez.
  def __apply_boundary_conditions(self):
    self.__Ez[0, :] = self.__Ez[-1, :] = 0
    self.__Ez[:, 0] = self.__Ez[:, -1] = 0
    # self.__Ez[~self.guide_mask] = 0


  def get_temp_diff(self) -> float:
    return self.__dt


  def set_grid_size(self, x: int, y: int):
    self.__grid_x = x
    self.__grid_y = y
    self.__Ez = np.zeros((self.__grid_x, self.__grid_y))  
    self.__Hx = np.zeros((self.__grid_x, self.__grid_y))  
    self.__Hy = np.zeros((self.__grid_x, self.__grid_y))


  def set_source(self, sources: list) -> None:
    self.__sources = sources


  def clear(self) -> None:
    self.__Ez.fill(0)
    self.__Hx.fill(0)
    self.__Hy.fill(0)
    self.__anim = None
    self.__ez_history.clear()


  def get_probe_values(self) -> list:
    return self.__ez_history


  def set_probe(self, x: int, y: int) -> None:
    self.__probe = (x,y)


  def run_no_graphics(self, n_steps) -> None:
    self.__n_steps = n_steps
    for t in range(self.__n_steps):
      self.__update_H()
      self.__update_E()
      self.__apply_boundary_conditions()
      self.__apply_source(t)
      self.__ez_history.append(self.__Ez[self.__probe[0], self.__probe[1]])


  def set_graphics(self, n_steps: int, frame_interval=30, mode=GraphicMode.Plane, mpl_cmap="inferno", interpol=None, title="FDTD Simulation") -> None:
    self.__n_steps = n_steps

    if mode == GraphicMode.Plane:
      fig, ax = plt.subplots(figsize=(8,6))
      im = ax.imshow(self.__Ez.T, origin='lower', cmap=mpl_cmap, vmin=-0.05, vmax=0.05, interpolation=interpol)
      ax.set_title(title)
      ax.set_xlabel("x (mm)")
      ax.set_ylabel("y (mm)")
      cbar = plt.colorbar(im, ax=ax, orientation='vertical')
      cbar.set_label('Campo elétrico Ez (V/m)')

      def animate(t):
        self.__update_H()
        self.__update_E()
        self.__apply_boundary_conditions()
        self.__apply_source(t)
        self.__ez_history.append(self.__Ez[self.__probe[0], self.__probe[1]])
        im.set_array(self.__Ez.T)
        ax.set_title(f"{title}\nTime: {(t * self.__dt * 2 * 1e9):.4f} ns")
        return [im]


    elif mode == GraphicMode.Mesh:
      fig = plt.figure(figsize=(10, 8))
      ax = fig.add_subplot(111, projection='3d')
      
      ax.set_box_aspect([
        self.__grid_x * self.__dx * 1e3,
        self.__grid_y * self.__dy * 1e3,
        self.__grid_y * self.__dy * 1e3
      ])

      X, Y = np.meshgrid(np.arange(self.__grid_x) * self.__dx * 1e3,
                          np.arange(self.__grid_y) * self.__dy * 1e3)

      # Holds reference to surface
      surf = [None]

      ax.set_zlim(-0.3, 0.3)
      ax.set_xlabel('x (mm)')
      ax.set_ylabel('y (mm)')
      ax.set_zlabel('Ez (V/m)')
      ax.set_title(title)

      def animate(t):
        self.__update_H()
        self.__update_E()
        self.__apply_boundary_conditions()
        self.__apply_source(t)
        self.__ez_history.append(self.__Ez[self.__probe[0], self.__probe[1]])

        # Removes the old surface to draw the new
        if surf[0] is not None:
            surf[0].remove()

        surf[0] = ax.plot_surface(X, Y, self.__Ez.T, cmap=mpl_cmap, edgecolor='none', vmin=-0.1, vmax=0.1)
        ax.set_title(f"{title}\nTime: {(t * self.__dt * 2 * 1e9):.4f} ns")
        return [surf[0]]

    self.__anim = animation.FuncAnimation(fig, animate, frames=self.__n_steps, interval=frame_interval, blit=True)



  def render(self, frame_rate=30, dpi=300, output_file="simulation.mp4") -> None:
    progress_bar = tqdm(total=self.__n_steps, desc="Rendering animation", unit="frame")    
    def progress_callback(i, n):
        progress_bar.update(1)
    self.__anim.save(output_file, writer="ffmpeg", fps=frame_rate, dpi=dpi, progress_callback=progress_callback)
    progress_bar.close()
