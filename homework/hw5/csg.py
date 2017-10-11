import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# =============================================================================
# =============================================================================
# Problem 1
# =============================================================================
# =============================================================================



# =============================================================================
# Point Class - Complete    
# =============================================================================

class Point :
    """A point is an (x, y) coordinate on the plane."""
    
    def __init__(self, x, y) :
        self.__x, self.__y = x, y
        
    @property
    def x(self) :
        return self.__x
    
    @property
    def y(self) :
        return self.__y
        
    def __str__(self) :
        """Prints the point"""
        return "(%.6f, %.6f)" % (self.__x, self.__y)
        
    def __add__(self, other) :
        """Adds two points to create a new point"""
        return Point(self.__x + other.x, self.__y + other.y)
        
    def __sub__(self, other) :
        """Adds two points to create a new point"""
        return Point(self.__x - other.x, self.__y - other.y)
        
    def __mul__(self, scale):
        """Scales a point by a scalar number"""
        new_x = self.__x * scale
        new_y = self.__y * scale
        return Point(new_x, new_y)
    
    def __truediv__(self, scale):
        """Scales a point by the inverse of a scalar number"""
        new_x = self.__x / scale
        new_y = self.__y / scale
        return Point(new_x, new_y)
    
    def distance(self, p = False) :
        if not p : 
            p = Point(0.0, 0.0)
        return np.sqrt((self.__x-p.x)**2 + (self.__y-p.y)**2)
    
# =============================================================================
# Ray Class - Given - Additions Made
# =============================================================================
      
class Ray:
    
    def __init__(self, origin, direction) :
        self.__origin = origin
        # ensure the direction is normalized to unity, i.e., cos^2 + sin^2 = 1
        norm = np.sqrt(direction.x**2 + direction.y**2)
        self.__direction = Point(direction.x/norm, direction.y/norm)
        self.__m = (origin.y - direction.y)/(origin.x - direction.x)
        self.__b = origin.y - self.__m * origin.x
            
    def __str__(self) :
        return "Ray: r_0(%10.6f, %10.6f), d(%.6f %.6f) " % \
               (self.origin.x, self.origin.y, self.direction.x, self.direction.y)
    @property
    def m(self) :
        return self.__m

    @property
    def b(self) :
        return self.__b   

    @property
    def origin(self) :
        return self.__origin
    
    @property
    def direction(self) :
        return self.__direction
               


# =============================================================================
# =============================================================================
# Problem 2
# =============================================================================
# =============================================================================



# =============================================================================
# Node Class - Abstract - Given
# =============================================================================
class Node:

    def contains(self, p) :
        """Does the node contain the point?"""
        raise NotImplementedError

    def intersections(self, r) :
        """Where does the node intersect the ray?"""
        raise NotImplementedError

# =============================================================================
# Primitave Class
# =============================================================================
class Primitive(Node) :
    
    def __init__(self, surface, sense) :
        """ Define a node consisting of a directed surface.
        Here, sense indicates "into" or "outof" the surface.
        
        Arguments:
            surface : Surface
                Surface (or derived variant) to define the node
            sense : bool
                True for inside the surface, False for outside 
                the surface
        """
        self.surface, self.sense = surface, sense
        
    def contains(self, p) :
        return (self.surface.f(p) < 0) == self.sense
        
    def intersections(self, r) :
        ints = self.surface.intersections(r)
        ints.sort(key=lambda p: p.x)
        return ints
        
# =============================================================================
# Operator Class
# =============================================================================

class Operator(Node) :
    
    def __init__(self, L, R) :
        self.L, self.R = L, R
        # some super checking algorithm

    def contains(self, p) :
        raise NotImplementedError

    def intersections(self, r) :
        # get intersections with left and right nodes
        pointsL = self.L.intersections(r)
        pointsR = self.R.intersections(r)
        ints = pointsL + pointsR
        ints.sort(key=lambda p: p.x)
        return ints
      
# =============================================================================
# Union Class
# =============================================================================
class Union(Operator):
    
    def __init__(self, L, R) :
        super(Union, self).__init__(L, R)
        
    def contains(self, p) :
        inL = self.L.contains(p)
        inR = self.R.contains(p)
        return inL or inR

# =============================================================================
# Intersection Class
# =============================================================================
class Intersection(Operator) :

    def __init__(self, L, R) :
        super(Intersection, self).__init__(L, R)

    def contains(self, p) :
        inL = self.L.contains(p)
        inR = self.R.contains(p)
        return inL and inR  


# =============================================================================
# =============================================================================
# Problem 3
# =============================================================================
# =============================================================================


      
# =============================================================================
# Surface Class        
# =============================================================================
class Surface :
    
    def f(self, p) :
        raise NotImplementedError
        
    def intersections(self, r) :
        raise NotImplementedError
        
# =============================================================================
# QuadraticSurface Class         
# =============================================================================
class QuadraticSurface(Surface) :
    
    def __init__(self, A=0.0, B=0.0, C=0.0, D=0.0, E=0.0, F=0.0) :
        """Surface is defined by the polynomial :
            f(x,y) = A*x**2 + B*y**2 + C*x*y + D*x + E*y + F 
            
        Arguments:
            A: float
                coefficient
            B: float
                coefficient
            C: float
                coefficient
            D: float
                coefficient
            E: float
                coefficient
            F: float
                coefficient
        """    
        self.__A = A
        self.__B = B
        self.__C = C
        self.__D = D
        self.__E = E
        self.__F = F
        self.__coeff = np.array([A, B, C, D, E, F])

    @property
    def A(self) :
        return self.__A

    @property
    def B(self) :
        return self.__B

    @property
    def C(self) :
        return self.__C

    @property
    def D(self) :
        return self.__D

    @property
    def E(self) :
        return self.__E

    @property
    def F(self) :
        return self.__F      
    
    @property
    def coeff(self) :
        return self.__coeff 
    
    def __str__(self) :
        return "%.3f*x**2 + %.3f*y**2 + %.3f*x*y + %.3f*x + %.3f*y + %.3f" % tuple(self.coeff)
    
    def intersections(self, r) :
        """Determines the first point where a ray intersects the geometry
        
        Arguments:
            self: object
                the QuadraticSurface the ray will intersect
            r: object of ray class
                the ray that will intersect the QuadraticSurface
            
        Returns:
            p: list 
                points of intersection
        """
        # Verify that r is a ray, otherwise raise TypeError
        try :
            assert type(r) == Ray
        except :
            raise TypeError("r must be a Ray")
        
        # Ray Eqn
        # y = mx + b
        # 0 = mx - y + b
        
        # System to sovle 
        # A*x**2 + B*y**2 + C*x*y + D*x + E*y + F = 0   (1)
        #                           m*x       + b = y   (2)
        #
        # Substitutin in Eq. 2 for y in Eq. 1 :
        # (A + B*m**2 + C*m)*x**2 + (2*B*b*m + C*b + D + E*m)*x + B*b**2 + E*b + F        (3)
        # If A = B = C = 0, the equation is linear :
        #   x = -(E*b + F)/(D + E*m)
        # Else the equation is quadratic :
        #   x1 = (-2*B*b*m - C*b - D - E*m + sqrt(-4*A*B*b**2 - 4*A*E*b - 4*A*F + 4*B*D*b*m - 4*B*F*m**2 + C**2*b**2 + 2*C*D*b - 2*C*E*b*m - 4*C*F*m + D**2 + 2*D*E*m + E**2*m**2))/(2*(A + B*m**2 + C*m))
        #   x2 = -(2*B*b*m + C*b + D + E*m + sqrt(-4*A*B*b**2 - 4*A*E*b - 4*A*F + 4*B*D*b*m - 4*B*F*m**2 + C**2*b**2 + 2*C*D*b - 2*C*E*b*m - 4*C*F*m + D**2 + 2*D*E*m + E**2*m**2))/(2*A + 2*B*m**2 + 2*C*m)
        
        A, B, C, D, E, F = self.coeff
        m, b = r.m, r.b
        
        # If the equation is not quadratic
        if not (A or B or C) :
            x = -(E*b + F)/(D + E*m)
            y = m*x + b
            return [Point(x, y)]
        # Else it is quadratic
        else :
            # If the sqrt contains a positive number, there are two intersections.
            if -4*A*B*b**2 - 4*A*E*b - 4*A*F + 4*B*D*b*m - 4*B*F*m**2 + C**2*b**2 + 2*C*D*b - 2*C*E*b*m - 4*C*F*m + D**2 + 2*D*E*m + E**2*m**2 > 0 :
                x1 = (-2*B*b*m - C*b - D - E*m + np.sqrt(-4*A*B*b**2 - 4*A*E*b - 4*A*F + 4*B*D*b*m - 4*B*F*m**2 + C**2*b**2 + 2*C*D*b - 2*C*E*b*m - 4*C*F*m + D**2 + 2*D*E*m + E**2*m**2))/(2*(A + B*m**2 + C*m))
                x2 = -(2*B*b*m + C*b + D + E*m + np.sqrt(-4*A*B*b**2 - 4*A*E*b - 4*A*F + 4*B*D*b*m - 4*B*F*m**2 + C**2*b**2 + 2*C*D*b - 2*C*E*b*m - 4*C*F*m + D**2 + 2*D*E*m + E**2*m**2))/(2*A + 2*B*m**2 + 2*C*m)
                y1 = m*x1 + b
                y2 = m*x2 + b
                ints = [Point(x1,y1), Point(x2,y2)]
                ints.sort(lambda p: p.x)
                return ints
            # Else if the sqrt contains a negative number, there are no intersections.
            elif -4*A*B*b**2 - 4*A*E*b - 4*A*F + 4*B*D*b*m - 4*B*F*m**2 + C**2*b**2 + 2*C*D*b - 2*C*E*b*m - 4*C*F*m + D**2 + 2*D*E*m + E**2*m**2 < 0 :
                return []
            # Else the sqrt contains zero, and there is one (tangent) intersection.
            else :
                x = (-2*B*b*m - C*b - D - E*m)/(2*(A + B*m**2 + C*m))
                y = m*x + b
                return [Point(x,y)]
                
    def f(self, p) :
        x, y = p.x, p.y
        A, B, C, D, E, F = self.coeff
        return A*x**2 + B*y**2 + C*x*y + D*x + E*y + F 


class PlaneV(QuadraticSurface) :
    
    def __init__(self, a) :
        """Subclass of QuatraticSurface where the surface is a 
        verical plane of the following equations :
            x = a
            
        Arguments:
            a: float
                x location of plane
        """
        super(PlaneV, self).__init__(D=1.0, F=-a)
        self.__a = a
        
    @property
    def a(self) :
        return self.__a


class PlaneH(QuadraticSurface) :
    
    def __init__(self, b) :
        """Subclass of QuatraticSurface where the surface is a 
        horizontal plane of the following equations :
            y = b
            
        Arguments:
            b: float
                y location of plane
        """
        super(PlaneH, self).__init__(E=1.0, F=-b)
        self.__b = b
        
    @property
    def b(self) :
        return self.__b
        

class Plane(QuadraticSurface) :
    
    def __init__(self, m, b) :
        """Subclass of QuatraticSurface where the surface is a 
        plane of the following equations :
            y = m*x + b
            
        Arguments:
            m: float
                slope of plane
            b: float
                y-intercept of plane
        """
        super(Plane, self).__init__(D=-m, E=1.0, F=-b)
        self.__m = m
        self.__b = b
        
    @property
    def m(self) :
        return self.__m
    
    @property
    def b(self) :
        return self.__b

class Circle(QuadraticSurface) :
    
    def __init__(self, r, a=0.0, b=0.0) :
        """Subclass of QuatraticSurface where the surface is a 
        circle of the following equations :
            x**2 + y**2 - 2*a*x - 2*b*y  + (-r**2 + b**2 + a**2) = 0
            
        Arguments:
            r: float
                radius of the circle
            a: float
                x coordinate of the circle's center
            b: float
                y coordinate of the circle's center
        """
        super(Circle, self).__init__(A=1, B=1, D=-2*a, E=-2*b, F=(-r**2 + b**2 + a**2))
        self.__r = r
        self.__a = a
        self.__b = b
        
    @property
    def r(self) :
        return self.__r
    
    @property
    def a(self) :
        return self.__a
    
    @property
    def b(self) :
        return self.__b
        
        
# =============================================================================
# =============================================================================
# Problem 4
# =============================================================================
# =============================================================================


               
# =============================================================================
# Region Class               
# =============================================================================
class Region :
    
    def __init__(self) :
        self.node = None
    
    def append(self, node=None, surface=None, operation="U", sense=False) :
        """Adds a node through Union or Intersection to the Region 
        based on a pased node or passed surface and sense.
        
        Arguments:
            
            node: Node
                node to be added 
            surface: Surface
                surface to made into node to be added
            operation: str
                "U" to perform Union. Other to perform Intersection
            sense: bool
                sense to be used to create node from surface
        """
        # Returns AssertionError if both node and surface are given
        # this means the function can only add a node or a surface
        # at a time.
        assert((node and not surface) or (surface and not node))
        # If surface is a surface
        if isinstance(surface, Surface) :
            # Defines a node based on surface with sense = sense
            node = Primitive(surface, sense)
        # If there is no existing node
        if self.node is None :
            # Define node to be node (either passed to the function or
            # based on the surace passed)
            self.node = node
        # Else, the new node must be added to the existing node
        else :
            # Defines O as a class Union object if operation="U"
            # Else defines O as class Intersection object.
            O = Union if operation == "U" else Intersection
            # Performs either the Union or Intersection of the self.node
            # and node
            self.node = O(self.node, node)
          
    def intersections(self, r) :
        return self.node.intersections(r)
    
    def contains(self, p) :
        return self.node.contains(p)
    

# =============================================================================
# Geometry Class        
# =============================================================================
class Geometry :
    
    # Attributes can be defined in the body of a class.  However, these
    # become "static" values that are the same for every object of the class.
    # Hence, they can be accessed either through object.attribute or 
    # classname.attribute.
    noregion = -1    
    
    def __init__(self,  xmin, xmax, ymin, ymax) :
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.regions = []
        
    def add_region(self, r) :
        self.regions.append(r)

    def find_region(self, p) :
        region = Geometry.noregion
        # look for the region containing p.
        return region
        
    def plot(self, nx, ny) :
        pass
        
if __name__ == '__main__' :
    
    # unit circle centered at origin
    c0 = Circle(1)
    # circle of radius two centered at the origin
    c1 = Circle(2)
    
    # produce a region that represents the area between the two circles
    region = Region()
    region.append(surface=c0, sense=False, operation="I")
    region.append(surface=c1, sense=True, operation="I")
    
    ray = Ray(Point(-3, 0), Point(1, 0))        
    ints = region.intersections(ray)

    
    
    
    