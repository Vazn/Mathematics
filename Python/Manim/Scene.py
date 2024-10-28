from manim import *
from sympy import *
import numpy as num

x, y, t = symbols("x, y, t")
F = Matrix([y, -sin(x - num.pi)])
f = lambdify([x, y, t], F)
print(*(f(1, 1, 0)[0]))
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
      background = NumberPlane()
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



