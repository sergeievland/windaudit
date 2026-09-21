import math
from dataclasses import replace
import pytest
from windaudit.frame import Node,NodePoint
from windaudit.calibrate import FrozenConfig
from windaudit.matching import MatchResult
from windaudit.order_constraints import build_constraints,negative_cycle,configuration


def config():
    return FrozenConfig(20.,18.,22.,20,1.,7.,14.,40.,math.radians(4),math.radians(120),2,'CW')

def point(pid,r,theta,label=0):
    return NodePoint(pid,[r*math.cos(theta),r*math.sin(theta),0.],r,theta,0.,label,0,label)

def empty():
    return MatchResult([],{}, {},[],0,0,0,0,0,0)

@pytest.mark.parametrize('gap,difference,expected',[(50.,3,2),(-50.,3,-4)])
def test_order_sign(gap,difference,expected):
    a=('relative',1);b=('relative',2)
    pairs=[(a,point(i,100+gap,1.,difference),b,point(i,100,1.),gap) for i in [1,2]]
    arcs,_=build_constraints({},config(),1.,empty(),pairs)
    assert len(arcs)==1 and arcs[0].c==expected
    assert (arcs[0].u,arcs[0].v)==((a,b) if gap>0 else (b,a))

def test_cut_transport_and_threshold():
    a=('relative',1);b=('relative',2)
    pairs=[(a,point(i,150,2*math.pi-.01),b,point(i,100,.01),50.) for i in [1,2]]
    arcs,_=build_constraints({},config(),1.,empty(),pairs)
    assert arcs[0].c==0 # d=+1, not the uncorrected d=0
    assert not build_constraints({},config(),2.5,empty(),pairs)[0]
    assert configuration(config(),1.)['sha256']!=configuration(config(),2.)['sha256']

def test_smooth_nonintersecting_sheet_does_not_imply_cross_angle_order():
    # r(t) = 1000 + D*t/(2*pi) + A*sin(20*t). Every successive
    # revolution lies exactly D farther out along the same ray; the sheet
    # never intersects itself. Within a four-degree sector, the SAME turn
    # nevertheless changes radius by more than D.
    D=20.;A=100.;t=1.2;dt=.05
    r=lambda x:1000+D*x/(2*math.pi)+A*math.sin(20*x)
    assert abs(r(t+dt)-r(t))>D
    for i in range(100):
        theta=i*2*math.pi/100
        assert r(theta+2*math.pi)-r(theta)==pytest.approx(D)
    a=('same',1);b=('same',2)
    gap=r(t)-r(t+dt)
    pairs=[(a,point(i,r(t),t),b,point(i,r(t+dt),t+dt),gap) for i in [1,2]]
    arcs,_=build_constraints({},config(),1.,empty(),pairs)
    assert len(arcs)==1
    # True same-turn gauges are both zero and violate the inferred bound.
    assert arcs[0].c==-1
