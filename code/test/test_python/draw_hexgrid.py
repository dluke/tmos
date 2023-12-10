
"""
pygame is not compatible with python3.8 (latest is python3.7, 2020-08)

As as a result this program is in an awkward spot. 
Perhaps we should be using python3.7 since it is compatible with MPL 3.3 and pygame
"""

import numpy as np
import pygame
pygame.init()


from tmos.base import Vector3d
import tmos.surface as surface

norm = np.linalg.norm

# colors
Color = pygame.Color
black = Color(0,0,0)
white = Color(255,255,255)
red = Color(255,0,0)
yellow = Color(255,255,0)
green = Color(0,255,0)
blue = Color(0,0,255)


scale = 100
ssize = (1024,1024)

offset = np.array(ssize)/2

########################################

# global mapping to and from the draw surface
def unmap(pt):
    return (pt - offset)/scale

def screenmap(pt):
    return scale * pt  + offset


# program flow
def _check_escape(events):
    touchedscape = lambda e: hasattr(e,'key') and e.key == pygame.K_ESCAPE
    return any(map(touchedscape, events))
def _check_quit(events):
    is_quit = pygame.QUIT in [e.type for e in events]
    is_escape = _check_escape(events)
    return is_quit or is_escape


#########################################
def hexc_to_array(hx):
    return np.array([hx[0],hx[1]])

def hexcseq_to_array(hlist):
    arr = np.zeros((len(hlist),2))
    for i, hx in enumerate(hlist):
        arr[i,0] = hx[0]
        arr[i,1] = hx[1]
    return arr

def v3d_to_npy2d(v):
    return np.array([v.x, v.y])

def npy2d_to_v3d(pt):
    x, y = pt
    return Vector3d(x, y, 0)

def make_xy_map(basis):
    def xy_map(hx):
        return hx[0]*basis[0] + hx[1]*basis[1]
    return xy_map


########################################
def construct(sphere_radius):
    hxsph = surface.HexSphereGrid(sphere_radius)
    return hxsph

def test_coordinates():
    R = 0.5
    hxsph = construct(R)
    vy = 2*R * Vector3d(0, 1, 0)
    hy = hxsph.get_rq(vy)
    vx = 2*R * Vector3d(np.sqrt(3.)/2., -1./2., 0)
    print(vx)
    print(vy) 
    hx = hxsph.get_rq(vx)
    print(hexc_to_array(hx))
    print(hexc_to_array(hy))

    print(hxsph.get_xyz(hx))
    print(hxsph.get_xyz(hy))

    # far away
    farvx = 100 * vx
    farhx = hxsph.get_rq(farvx)
    print(hexc_to_array(farhx))


#########################################

def _new_surface():
    screen = pygame.display.set_mode(ssize, 0, 32)
    return screen 

class Grabbable(object):

    color = blue
    gcolor = yellow
    radius = 5

    def __init__(self, pt):
        self.pt = pt
        self.grabbed = False

    def grab(self):
        self.grabbed = True

    def release(self):
        self.grabbed = False

    def update(self, rel):
        if self.grabbed:
            self.pt += rel 

    def draw(self, screen):
        color = self.gcolor if self.grabbed else self.color
        pygame.draw.circle(screen, color, np.rint(self.pt).astype(int), self.radius, 0)

    def _mouse_over(self, pos):
        return norm(np.array(pos) - self.pt) < self.radius

    def take_input(self, events):
        for ev in events:
            if ev.type == pygame.MOUSEBUTTONDOWN and self._mouse_over(ev.pos):
                # clicked on this grabbable
                self.grab()
            if ev.type == pygame.MOUSEBUTTONUP and self._mouse_over(ev.pos):
                self.release()
            if ev.type == pygame.MOUSEMOTION and self.grabbed:
                self.update(ev.rel)

class Body(object):

    color = white

    def __init__(self, head, tail, R, length=None):
        self.head = Grabbable(head)
        self.tail = Grabbable(tail)
        self.R = R
        self.length = length

    def update_length(self):
        self.length = np.linalg.norm(self.head.pt -self.tail.pt)

    def get_axis(self):
        return npy2d_to_v3d(unmap(self.head.pt) - unmap(self.tail.pt)).unit()

    def get_center(self):
        return (unmap(self.head.pt) + unmap(self.tail.pt))/2

    def get_capsule(self):
        self.update_length()
        center = self.get_center()
        caps = surface.Capsule(npy2d_to_v3d(center),
                self.get_axis(),
                self.R,
                self.length/scale)
        print(caps)
        return caps

    def draw(self, screen):
        self.head.draw(screen)
        self.tail.draw(screen)
        width = 2
        R = int(self.R * scale)
        pygame.draw.circle(screen, self.color, np.rint(self.head.pt).astype(int), R, width)
        pygame.draw.circle(screen, self.color, np.rint(self.tail.pt).astype(int),  R, width)
        midpt = (self.head.pt + self.tail.pt)/2
        pygame.draw.circle(screen, green, np.rint(midpt).astype(int), 5, 0)

    def is_grabbed(self):
        return self.head.grabbed or self.tail.grabbed


    def take_input(self, events):
        self.head.take_input(events)
        self.tail.take_input(events)


r3 = np.sqrt(3)
R60 = np.array([[1./2, -r3/2.], [r3/2., 1./2]])

def construct_polygon():
    start = np.array([1./(2*r3), 1./2])
    move = np.array([-1./r3, 0])
    poly = np.zeros((6,2))
    poly[0] = start
    for i in range(1,6):
        step = poly[i-1] + move
        poly[i] = step
        move = np.dot(R60, move)
    return poly

class State(object):

    def __init__(self):
        self.poly = construct_polygon()
        #self.disp = pygame.display.get_surface.get_size()
        self.body = None
        self.hxsph = None

    def update(self, events):
        # hangle events
        self.body.take_input(events)
        if not self.body.is_grabbed():
            return 
        self.update_body()

    def update_body(self):
        # get highlighting
        caps = self.body.get_capsule()
        # returns hexc coordinates 
        hxs = self.hxsph.body_set(caps)
        hxs = list(map(self.hxsph.get_xyz, hxs))
        hxy = list(map(v3d_to_npy2d, hxs))
        self.drawables['hl'] = self._highlight(hxy)

    def _highlight(self, hxy):
        hxscale = scale * self.hxsph.get_R/0.5
        h_hexes = [(hxscale*self.poly) + screenmap(pt) for pt in hxy]
        return h_hexes


def draw_pts(screen, xy, style):
    radius = style.get('radius')
    color = style.get('color')
    disp = np.array(screen.get_size())/2
    for pt in xy:
        spt = (scale * pt).astype(int) + disp
        pygame.draw.circle(screen, color, spt, radius, 0)

def draw_hexes(screen, hexes, style):
    color = style.get('color')
    width = style.get('width')
    disp = np.array(screen.get_size())/2
    for poly in hexes:
        pygame.draw.polygon(screen, color, poly, width)
        
def draw(screen, drawables):
    hlstyle = {'color':red, 'width':0}
    if 'hl' in drawables:
        draw_hexes(screen, drawables['hl'], hlstyle)
    ptstyle = {'radius':2, 'color':white}
    hxstyle = {'color':white, 'width':1}
    draw_pts(screen, drawables['xy'], ptstyle)
    draw_hexes(screen, drawables['hexes'], hxstyle)
    drawables['body'].draw(screen)


def main():
    sphere_radius = 0.25
    hxsph = construct(sphere_radius)
    poly =  construct_polygon()

    origin = surface.Hexc(0,0)
    # returns list of Hexc objects, range = 
    size = 6
    hlist = hxsph.coordinate_range(origin, size)
    harr = hexcseq_to_array(hlist)
    basis = list(map(v3d_to_npy2d, (hxsph.get_hrr(), hxsph.get_hrq())))
    xy_map = make_xy_map(basis)

    screen = _new_surface()
    xy = list(map(xy_map, harr))

    hxscale = scale * sphere_radius/0.5
    hexes = [(hxscale*poly) + screenmap(pt) for pt in xy]
    
    import collections
    drawables = collections.OrderedDict()
    drawables['xy'] = xy
    drawables['hexes'] = hexes
    disp = np.array(screen.get_size())/2
    tmphead = screenmap(np.array([1.,0.]))
    tmptail = screenmap(np.array([-1.,0.]))
    #print 'mapping (1,0) to', tmphead
    R = 0.5

    body = Body(tmphead, tmptail, 0.5)
    drawables['body'] = body
    state = State()
    state.hxsph = hxsph
    state.body = body
    state.drawables = drawables
    state.update_body()

    clock = pygame.time.Clock()
    fps = 40
    while True:
        # event handling
        screen.fill(black)
        events = pygame.event.get()
        if _check_quit(events): break
        if len(events) != 0:
            state.update(events)
        # draw
        draw(screen, drawables)
        pygame.display.flip()

        clock.tick(fps)

if __name__=='__main__':

    main()
    test_coordinates()
