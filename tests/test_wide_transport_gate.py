"""The bulk parent lookup must preserve the original integrability audit."""
import numpy as np
import pytest
from windaudit.scanspace import ScanPatch, audit_edge_integrability
from windaudit.wide_transport_gate import audit

def origin(z):return np.zeros(np.shape(z)+(2,))

@pytest.mark.parametrize('kind',['plane','hole','disconnected'])
def test_bulk_gate_matches_reference(kind):
    y,x=np.mgrid[-5:6,-5:6]
    xyz=np.stack([x+.23,y+.17,np.ones_like(x)],axis=-1)
    if kind=='plane':xyz[...,0]+=20
    elif kind=='hole':xyz[4:7,4:7]=-1
    else:xyz[:,4:7]=-1
    patch=ScanPatch(xyz)
    assert audit(patch,origin)==audit_edge_integrability(patch,origin)
