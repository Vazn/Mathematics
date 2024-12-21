from manim import *
from sympy import *
import numpy as num

AXISCONFIG = {
    'x_range':[-25, 25, 1],
    'y_range':[-25, 25, 1],
    'x_length':50,
    'y_length':50,
    "color": '#E8E8E8',
    "stroke_width": 4,
    "stroke_opacity": 1,
}
PLANECONFIG = {
    'x_range':[-25, 25, 1],
    'y_range':[-25, 25, 1],
    'x_length':50,
    'y_length':50,
    "color": '#E8E8E8',
    "stroke_width": 4,
    "stroke_opacity": 1,

    'faded_line_ratio' :2,
    'background_line_style':{
        "stroke_width": 3,
        "stroke_opacity": 0.35,
        "stroke_color": '#FFFFFF'
    },
    'faded_line_style':{
        "stroke_width": 3,
        "stroke_opacity": 0.10,
        "stroke_color": '#FFFFFF'
    },
}
POLARPLANECONFIG = {
    'size': 30,
    'radius_max': 15,
    'radius_step':1,
    'azimuth_step': 8,
    'faded_line_ratio' :3,
    'background_line_style':{
        "stroke_width": 4,
        "stroke_opacity": 0.95,
        "stroke_color": '#6E8CB6'
    },
    'faded_line_style':{
        "stroke_width": 3,
        "stroke_opacity": 0.3,
        "stroke_color": '#6E8CB6'
    },
}
SPACECONFIG = {}

x, y, z, t = symbols("x, y, z, t")
F = Matrix([y, -sin(x - num.pi)])
f = lambdify([x, y, t], F)
time = ValueTracker(0)

class VectorFlow(Scene):
   def createParticles(self, density):
      coeff = 1/density
      Q1 = [Dot(radius= coeff/3, point=[0 + i*coeff, 0 + j*coeff, 0]) for i in range(density) for j in range(density)]
      Q2 = [Dot(radius= coeff/3, point=[0 + i*coeff, 0 - j*coeff, 0]) for i in range(density) for j in range(density)]
      Q3 = [Dot(radius= coeff/3, point=[0 - i*coeff, 0 - j*coeff, 0]) for i in range(density) for j in range(density)]
      Q4 = [Dot(radius= coeff/3, point=[0 - i*coeff, 0 + j*coeff, 0]) for i in range(density) for j in range(density)]

      return list(set(Q1 + Q2 + Q3 + Q4))
   
   def construct(self):
      background = NumberPlane(**PLANECONFIG)
      vectorField = always_redraw(
         lambda : ArrowVectorField(
            lambda pos: f(pos[0], pos[1], time.get_value())[0] * RIGHT + f(pos[0], pos[1], time.get_value())[1] * UP
         )
      )
      particles = self.createParticles(20)
      
      self.add(background, vectorField)
      self.wait(1)

      self.add(*particles)
      self.wait(1)

      for particle in particles: particle.add_updater(vectorField.get_nudge_updater())
      self.play(time.animate.set_value(2), rate_func=linear, run_time=2)

class Other(Scene):
   pass
