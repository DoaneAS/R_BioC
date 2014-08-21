from swampy.TurtleWorld import *
import math
world = TurtleWorld()
bob = Turtle()
def square(t):
  for i in range(4):
    fd(t,100)
    lt(t)
def polygon(t, n, length):
    angle = 360.0 / n
    for i in range(n):
        fd(t, length)
        lt(t, angle)
def circle(t, r):
    circumference = 2 * math.pi * r
    n = int(circumference / 3) + 1
    length = circumference / n
    polygon(t, n, length)
print bob
circle(bob, 10)
wait_for_user()