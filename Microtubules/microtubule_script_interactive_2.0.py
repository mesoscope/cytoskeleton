# -*- coding: utf-8 -*-
import sys,os
import upy,numpy
upy.setUIClass()

from upy import uiadaptor

helperClass = upy.getHelperClass()
helper = helperClass()

class microtubuleUI(uiadaptor):
    def setup(self,vi=None):
        #self.helper = helperClass(vi=vi)
        self.dock = False
        self.spline = None
        self.initWidget()
        self.setupLayout()
        defaultpts = ((0,0,0), (0,0,500), (50,0,1000), (0,50,1500))
        self.points = defaultpts
        self.microtubule_model(defaultpts)
        self.synchronized = False
        self.current = None

    def CreateLayout(self):
        self._createLayout()
        return 1

    def Command(self, *args):
        self._command(args)
        return 1

    def initWidget(self, id=None):
        self.labelSel = self._addElemt(label = "Please select a spline to use for the model",width = 80,height = 10)
        self.selectedSpline = self._addElemt(name = "selSpline", width = 80, height = 10, action = None, type = "inputStr", variable = self.addVariable("str",''))
        self.synchronize = self._addElemt(name = "realTime", width=80,height=10,
                                             action=self.setSynchro,type="checkbox",
                                             variable=self.addVariable("int",0),
                                             label="Update in real time (slow! not recommended)")
        self.mtmake = self._addElemt(name="make",width=80,height=10,
                                     action=self.choose_spline,type="button",icon=None,
                                     variable=self.addVariable("int",0),label="Go")

    def setupLayout(self, ftype = "frame"):
        self._layout = []
        self._layout.append([self.labelSel])
        self._layout.append([self.selectedSpline])
        self._layout.append([self.mtmake])
        self._layout.append([self.synchronize])


    def choose_spline(self):
        helper.deleteChildrens(helper.getObject("microtubule"))
        self.current = self.getVal(self.selectedSpline)
        self.points = helper.getMeshVertices(helper.getObject(self.current),selected = False,transform = True)
        #self.points = helper.DecomposeMesh(helper.getObject(self.current))[1]
        print helper.getObject(self.getVal(self.selectedSpline)), 'choose'
        self.microtubule_model(self.points)
    def setSynchro(self,*args):
        if not self.synchronized:
            helper.synchronize(self.update)
            self.synchronized = True
        elif self.synchronized: 
            helper.unsynchronize(self.update)
            self.synchronized = False

    def update(self,*args,**kw):
        if self.selectedSpline is not None:
            print helper.getObject(self.getVal(self.selectedSpline)), 'update'
            checkpts = helper.getMeshVertices(helper.getObject(self.current), selected = False,transform = True)
            #checkpts = helper.DecomposeMesh(helper.getObject(self.current))
            if checkpts != self.points:
                helper.deleteChildrens(helper.getObject("microtubule"))
                self.microtubule_model(checkpts)
                self.points = checkpts
    #
    # Make a geometric model of a microtubule in upy
    #
    def microtubule_model(self,path,                             # points used for cubic spline
                          radius = 90.0, 	         	# microtubule center to tubulin center (Angstroms)
                          monomer_spacing = 40.0,	  	# distance along protofilament (Angstroms)
                          number_of_protofilaments = 13,	
                          rise_per_turn = -3,		# number of monomers (negative for left hand helix)
                          colors = ((.7,.7,.6,1),(.4,.4,.5,1))):

      spt = smooth_path(path, monomer_spacing)	  	# points and tangents
      pt = equispaced_points(spt, monomer_spacing)
      f = path_coordinate_frames(pt)
      pl = helix_subunit_placements(f, radius, monomer_spacing, number_of_protofilaments, rise_per_turn)
      r = 0.5*monomer_spacing
      col = [colors[0]]*number_of_protofilaments + [colors[1]]*number_of_protofilaments
      s = place_spheres(pl, r, col)
      #s.name = 'microtubule'
      return s

def smooth_path(path, spacing):

  #from Matrix import distance
  seg_lengths = [helper.measure_distance(path[i],path[i+1]) for i in range(len(path)-1)]
  subdiv = max(0, int(max(seg_lengths) / spacing) - 1)
  #from VolumePath import spline
  pt = natural_cubic_spline(path, subdiv, return_tangents = True)
  return pt

def natural_cubic_spline(points, segment_subdivisions, return_tangents = False):

  n = len(points)
  if isinstance(segment_subdivisions, int) and n > 0:
    segment_subdivisions = [segment_subdivisions] * (n - 1)

  d = len(points[0]) if n > 0 else 0

  # TODO: use tridiagonal solver to save time/memory.  SciPy.linalg.solve_banded
  from numpy import zeros, float32, empty, linalg
  a = zeros((n,n), float32)
  b = zeros((n,), float32)
  for i in range(1,n-1):
    a[i,i] = 4
    a[i,i-1] = a[i,i+1] = 1
  a[0,0] = a[n-1,n-1] = 1

  c = n + sum(segment_subdivisions)
  p = empty((c,d), float32)
  if return_tangents:
    tan = empty((c,d), float32)
  if n == 1:
    p[0,:] = points[0]
    if return_tangents:
      tan[0,:] = 0
  else:
    for axis in range(d):
      for i in range(1,n-1):
        b[i] = points[i+1][axis] -2*points[i][axis] + points[i-1][axis]
      z = linalg.solve(a, b)
      k = 0
      for i in range(n-1):
        div = segment_subdivisions[i]
        pc = div + 1 if i < n-2 else div + 2
        for s in range(pc):
          t = s / (div + 1.0)
          ct = points[i+1][axis] - z[i+1]
          c1t = points[i][axis] - z[i]
          p[k,axis] = z[i+1]*t**3 + z[i]*(1-t)**3 + ct*t + c1t*(1-t)
          if return_tangents:
            tan[k,axis] = 3*z[i+1]*t**2 - 3*z[i]*(1-t)**2 + ct - c1t
          k += 1

  if return_tangents:
    #import Matrix
    tan = helper.unit_vector(tan, axis = 1)
    return zip(p, tan)

  return p


    
def equispaced_points(points_and_tangents, spacing):

  pt = points_and_tangents
  ept = [pt[0]]
#  from VolumePath import spline
  arcs = arc_lengths([p for p,t in pt])
  d = spacing
  #from Matrix import linear_combination, normalize_vector
  for i, a in enumerate(arcs):
    while a > d:
      f = (d - arcs[i-1]) / (arcs[i] - arcs[i-1])
      (p0,t0),(p1,t1) = pt[i-1],pt[i]
      p = linear_combination((1-f), p0, f, p1)
      t = helper.unit_vector(linear_combination((1-f), t0, f, t1), axis = 1)
      ept.append((p,t))
      d += spacing
  return ept

# -----------------------------------------------------------------------------
# Return a list of arc lengths for a piecewise linear curve specified by
# points.  The points should be objects with operator - defined such as
# NumPy arrays.  There number of arc lengths returned equals the
# number of points, the first arc length being 0.
#
def arc_lengths(points):

  import math
  arcs = [0]
  for s in range(len(points)-1):
    d = points[s+1] - points[s]
    length = math.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])
    arcs.append(arcs[s] + length)
  return arcs

# -----------------------------------------------------------------------------
# Return vector a*u + b*v where a, b are scalars.
#
def linear_combination(*fv):

  flist = fv[::2]
  vlist = fv[1::2]
  fvlist = zip(flist,vlist)
  n = len(vlist[0])
  vsum = tuple([sum([f*v[a] for f,v in fvlist]) for a in range(n)])
  return vsum


def path_coordinate_frames(pt):

  #from VolumePath import tube
  f = extrusion_transforms(pt)
  return f

# -----------------------------------------------------------------------------
# Compute transforms mapping (0,0,0) origin to points along path with z axis
# along path tangent.
#
def extrusion_transforms(ptlist, yaxis = None):

    #import Matrix as M
    tflist = []
    if yaxis is None:
        # Make xy planes for coordinate frames at each path point not rotate
        # from one point to next.
        tf = numpy.array(identity_matrix()) #multiply_matrices assumes numpy
        n0 = (0,0,1)
        for p1,n1 in ptlist:
            tf = multiply_matrices(vector_rotation_transform(n0,n1), tf)
            tflist.append(multiply_matrices(translation_matrix(p1),tf))
            n0 = n1
    else:
        # Make y-axis of coordinate frames at each point align with yaxis.
        for p,t in ptlist:
            za = t
            xa = helper.unit_vector(numpy.cross(yaxis, za), axis = 1)
            ya = numpy.cross(za, xa)
            tf = ((xa[0], ya[0], za[0], p[0]),
                  (xa[1], ya[1], za[1], p[1]),
                  (xa[2], ya[2], za[2], p[2]))
            tflist.append(tf)
    return tflist


def identity_matrix():
    return ((1.0,0,0,0), (0,1.0,0,0), (0,0,1.0,0))

# -----------------------------------------------------------------------------
#
def multiply_matrices(*tf_list):

  if len(tf_list) == 2:
    tf1, tf2 = tf_list
    r1 = tf1[:,:3]
    t1 = tf1[:,3]
    r2 = tf2[:,:3]
    t2 = tf2[:,3]
    from numpy import zeros, float, add
    from numpy import dot as matrix_multiply
    tf = zeros((3,4), float)
    r = tf[:,:3]
    t = tf[:,3]
    r[:] = matrix_multiply(r1, r2)
    t[:] = add(t1, matrix_multiply(r1, t2))
  else:
    tf = multiply_matrices(*tf_list[1:])
    tf = multiply_matrices(tf_list[0], tf)
  return tf

# -----------------------------------------------------------------------------
#
def vector_rotation_transform(u,v): # I think there are ways to replace this with upy
                                    # but at this point it's easier and better to stick
                                    # to chimera

  cuv = numpy.cross(u,v)
  axis = helper.unit_vector(cuv, axis = 1)
  from math import atan2, pi
  angle = atan2(helper.vector_norm(cuv),numpy.dot(u,v))
  r = helper.rotation_matrix(angle, axis, (0,0,0))
  rotation =  numpy.array((tuple(r[0]), tuple(r[1]), tuple(r[2])))
  return rotation
# -----------------------------------------------------------------------------
#
def translation_matrix(shift): #this might also be done in pure upy in some way or another

  from numpy import array
  tf = array(((1.0, 0, 0, shift[0]),
              (0, 1.0, 0, shift[1]),
              (0, 0, 1.0, shift[2])))
  return tf


def helix_subunit_placements(f, radius, subunit_spacing, subunits_per_turn, rise_per_turn):

  pl = []
  #import Matrix as M
  for tf in f:
    for i in range(subunits_per_turn):
      a = i * (360.0/subunits_per_turn) # degrees
      dz = subunit_spacing*i*float(rise_per_turn)/subunits_per_turn
      from math import pi
      r = helper.rotation_matrix((a * pi / 180), (0,0,1), (0,0,0))
      rot = numpy.array((tuple(r[0]), tuple(r[1]), tuple(r[2])))
      stf = multiply_matrices(tf,
                                translation_matrix((0,0,dz)),
                                rot,
                                translation_matrix((radius,0,0)))
      pl.append(stf)

  return pl


def place_spheres(pl, radius, colors = [(.7,.7,.7,1)]):
    mt = helper.newEmpty("microtubule")
    spherecount = 0
    for num, matrix in enumerate(pl):
        pos = chimeraToUpyMat(matrix)[3]

        #helper.Sphere(name = "sphere"+str(spherecount), radius = radius, res = 25, color = colors[num%len(colors)], pos = pos)
        #helper.setTransformation("sphere"+str(spherecount), chimeraToUpyMat(matrix))
        mat = helper.addMaterial(name = "cube" + str(spherecount), color = colors[num%len(colors)]) #the way cubes work in upy, functions like colormaterial won't help:
                                                                                                    #they will just make all the cubes the same color
        helper.box(name= "cube"+str(spherecount), center = pos, size = [40,40,40], mat = mat, parent = mt) #mt from cubes
        helper.setTransformation("cube"+str(spherecount), chimeraToUpyMat(matrix))
##        if spherecount == 202:
##            helper.box(name = 'fakemotor', size = [30,30,30],  parent= helper.getObject("cube202"))
##            helper.setTransformation('fakemotor',chimeraToUpyMat(matrix))
##
##            helper.setTranslation('fakemotor',[35,0,0], absolue = False) #we're translating along the z-axis in c4d's world.. x axis in our world
    
        spherecount = spherecount + 1

    return

def chimeraToUpyMat(m):
    #chimera uses 3x4 matrix, upy uses 4x4
    M=[[1.0,0.,0.0,0.0],
       [0.0,1.,0.0,0.0],
       [0.0,0.,1.0,0.0],
       [0.0,0.,0.0,1.0]]

    xaxis = []
    yaxis = []
    zaxis = []
    for i in range(len(m)):
        xaxis.append(m[i][0])
        yaxis.append(m[i][1])
        zaxis.append(m[i][2])

    M[0][:3] = helper.ToVec(xaxis)
    M[1][:3] = helper.ToVec(yaxis)
    M[2][:3] = helper.ToVec(zaxis)
    posVec = helper.ToVec([m[0][3],m[1][3],m[2][3]])
    M[3][:3] = posVec
    return numpy.array(M) #so that we can use upy's FromMat



if uiadaptor.host == "tk":
    from DejaVu import Viewer
    vi = Viewer()    
    #require a master
    #import Tkinter #Tkinter if python2.x tkinter for python3.x
    #root = Tkinter.Tk()
    mygui = microtubuleUI(title="microtubuleUI",master=vi)
    mygui.setup(vi=vi)
    #mygui.display()
elif uiadaptor.host == "qt":
    from PyQt4 import QtGui
    app = QtGui.QApplication(sys.argv)
    mygui = microtubuleUI(title="microtubule")
    mygui.setup()
    #ex.show()
else :
    mygui = microtubuleUI(title="microtubuleUI")
    mygui.setup()
    #call it
mygui.display()
if uiadaptor.host == "qt": app.exec_()#work without it ?
