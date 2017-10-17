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
        """If scalar is given, scales a point by the inverse of a scalar number.
        If a point is given, the x and y terms are divided.
        """
        if isinstance(scale, Point):
            new_x = self.__x / scale.x
            new_y = self.__y / scale.y
        else :
            try : 
                assert (type(scale) == float) or (type(scale) == int)
            except :
                raise TypeError('scale must be a Point, float, or int')
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
        """Creates a ray objects
        
        Arguments:
            origin: Point
                The origin point of the ray
            direction: Point
                A vector originating at the origin defining the direction
                or the ray.
        """
        self.__origin = origin
        # ensure the direction is normalized to unity, i.e., cos^2 + sin^2 = 1
        norm = np.sqrt(direction.x**2 + direction.y**2)
        try :
            assert norm != 0
        except :
            raise ValueError('direction must have magnitude')
        self.__direction = Point(direction.x/norm, direction.y/norm)
        self.__m = (direction.y)/(direction.x)
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
    
    def intersects(self, p) :
        """Determines if a ray intersects a point
        
        Arguments:
            p: Point
                The point for which intersection is determined
                
        Returns:
            intersect: bool
                True if intersected, Else false
        """
        
        # Move reference origin to origin of ray
        new_p = p - self.origin
        # If the ray has no magnitude in the x-direction
        if not self.direction.x :
            #If the adjusted point's x-value is also 0
            if self.direction.x == new_p.x :
                # get the magnitude from y points, if positive point is on ray
                return new_p.y / self.direction.y >= 0
            else :
                return False
        # If the ray has no magnitude in the y-direction
        elif not self.direction.y :
            #If the adjusted point's y-value is also 0
            if self.direction.y == new_p.y :
                # get the magnitude from x points, if positive point is on ray
                return new_p.x / self.direction.x >= 0
            else :
                return False
        # ray has magnitude in x and y direction
        else :
            # divide
            d_point = new_p / self.direction
            # magnitudes must be equal for point to be on ray
            if d_point.x == d_point.y :
                d = d_point.x
                # magnitude must be positive to be on ray
                return d >= 0
            else :
                return False
         
               


# =============================================================================
# =============================================================================
# Problem 2
# =============================================================================
# =============================================================================



# =============================================================================
# Node Class - Abstract - Given
# =============================================================================
class Node:
    """Abstract class outlining methods neccessary for sub classes
    
    subclasses:
        Primative
        Operator
    """

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
        """Defines a node consisting of a directed surface.
        Here, sense indicates inside or outside the surface.
        
        Super:
            Node
        
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
        ints.sort(key=lambda p: p.distance(r.origin))
        return ints
        
# =============================================================================
# Operator Class
# =============================================================================

class Operator(Node) :
    
    def __init__(self, L, R) :
        self.L, self.R = L, R
        try :
            assert isinstance(L, Node)
            assert isinstance(R, Node)
        except :
            raise TypeError('L and R must be instances of Node')

    def contains(self, p) :
        raise NotImplementedError

    def intersections(self, r) :
        # get intersections with left and right nodes
        pointsL = self.L.intersections(r)
        pointsR = self.R.intersections(r)
        ints = pointsL + pointsR
        ints.sort(key=lambda p: p.distance(r.origin))
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
            # If the slope of the ray and line are equal
            if -D/E == m :
                # They do not intersect
                return []
            else :
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
                ints = [i for i in ints if r.intersects(i)]
                ints.sort(key=lambda p: p.distance(r.origin))
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
        """Retruns value of function that difines the surface evaluated at
        point p.
        
        Arguments:
            p: Point
                point to evaluate funtion at
        
        Returns:
            f(p): float
                function evaluated at p
        """
        try :
            assert type(p) == Point
        except :
            raise TypeError('p must be a point')
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
        try :
            assert type(r) == Ray
        except :
            raise TypeError('r must be a ray')
        ints = list(set(self.node.intersections(r)))
        ints.sort(key=lambda p: p.distance(r.origin))
        return ints
    
    def contains(self, p) :
        try :
            assert type(p) == Point
        except :
            raise TypeError('p must be a point')
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
        self.__xmin, self.__xmax, self.__ymin, self.__ymax = xmin, xmax, ymin, ymax
        self.__regions = []
        
    @property
    def xmin(self) :
        return self.__xmin
    
    @property
    def xmax(self) :
        return self.__xmax
    
    @property
    def ymin(self) :
        return self.__ymin
    
    @property 
    def ymax(self) :
        return self.__ymax
    
    @property
    def regions(self) :
        return self.__regions
        
    def add_region(self, r) :
        self.regions.append(r)

    def find_region(self, p) :
        """Finds the index(es) of the region(s) containing the point given
        
        Arguments:
            p: Point
                Find regions that contain this point
                
        Returns:
            region: list or int
                list of the index(es) of the region(s) containing p. If no 
                region contains p or p is out of the bounds of the geometry,
                Geometry.noregion is returned.
        """
        try :
            assert type(p) == Point
        except :
            raise TypeError('p must be a point')
        # Checks if point is within the bounds of the geometry.
        if p.x < self.__xmin or p.x > self.__xmax or p.y < self.__ymin or p.y > self.__ymax :
            return Geometry.noregion
        region = []
        # Checks each region, if it contains the point, adds the index to 
        # region
        for i in range(len(self.__regions)) :
            if self.__regions[i].contains(p) :
                region.append(i)
        # if no regions contained the point, returns Geometry.noregion
        if not len(region) :
            region = Geometry.noregion
        return region
        
    def plot(self, nx, ny) :
        """Plots the geometry
        
        Arguments:
            nx: int
                number of x points
            ny: int
                number of y points
        """
        x = np.linspace(self.xmin, self.xmax, nx)
        y = np.linspace(self.ymin, self.ymax, ny)
        
        Z = np.zeros((nx, ny))
        
        
        for i in range(len(x)) :
            for j in range(len(y)) :
                regions = self.find_region(Point(x[i], y[j]))
                if type(regions) == list :
                    regions = len(regions)                
                Z[i][j] = regions
        levels = list(range(-1, int(np.amax(Z))+1))        
        cp = plt.contourf(x,y,Z,levels)
        plt.colorbar(cp)
        plt.title('Geometry Plot')
        plt.show()
    
        
if __name__ == '__main__' :
    pass

    

    
    
    
    