#!/bin/python

import numpy as np
from pysundials.core import *
from pysundials.kinsol import *
from problems import AnalyticNonlinearSys


def test_kinsol():
    NEQ = 3
    m_aa = 2
    tol = 1e-4

    sunctx = SUNContextView.Create()
    kin_view = KINView.Create(KINCreate(sunctx.get()))

    problem = AnalyticNonlinearSys()

    u = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))

    def fp_function(u, g, _):
        return problem.fixed_point_fn(u, g)

    kin_status = KINSetMAA(kin_view.get(), m_aa)

    kin_status = KINInit(kin_view.get(), fp_function, u.get())

    kin_status = KINSetFuncNormTol(kin_view.get(), tol)

    # initial guess
    udata = N_VGetArrayPointer(u.get())
    udata[:] = [0.1, 0.1, -0.1]

    # no scaling used
    scale = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))
    N_VConst(1.0, scale.get())

    kin_status = KINSol(kin_view.get(), u.get(), KIN_FP, scale.get(), scale.get())
    print(f"kin_status={kin_status}, solution={udata}")
    
    ucor = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))
    ucor_data = N_VGetArrayPointer(ucor.get())
    
    problem.solution(ucor.get())

    assert(np.allclose(udata, ucor_data))
    
    kin_status, nni = KINGetNumNonlinSolvIters(kin_view.get(), 0)
    kin_status, nfe = KINGetNumFuncEvals(kin_view.get(), 0)
    print("Final Statistics:")
    print(f"nni={nni}, nfe={nfe}")


if __name__ == "__main__":
    test_kinsol()
