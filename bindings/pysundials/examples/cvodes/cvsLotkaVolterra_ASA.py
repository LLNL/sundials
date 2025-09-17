import numpy as np
from pysundials.core import *
from pysundials.cvodes import *

NEQ = 2
NP = 4
T0 = 0.0
TF = 10.0
RTOL = 1e-10
ATOL = 1e-14
STEPS = 5

params = np.array([1.5, 1.0, 3.0, 1.0], dtype=np.float64)


# Lotka-Volterra ODE system
def lotka_volterra(t, y, ydot, p):
    y = N_VGetArrayPointer(y)
    ydot = N_VGetArrayPointer(ydot)
    ydot[0] = p[0] * y[0] - p[1] * y[0] * y[1]
    ydot[1] = -p[2] * y[1] + p[3] * y[0] * y[1]
    return 0


def vjp(v, Jv, t, y, p):
    v = N_VGetArrayPointer(v)
    Jv = N_VGetArrayPointer(Jv)
    y = N_VGetArrayPointer(y)
    Jv[0] = (p[0] - p[1] * y[1]) * v[0] + p[3] * y[1] * v[1]
    Jv[1] = -p[1] * y[0] * v[0] + (-p[2] + p[3] * y[0]) * v[1]
    return 0


def parameter_vjp(v, Jv, t, y, p):
    v = N_VGetArrayPointer(v)
    Jv = N_VGetArrayPointer(Jv)
    y = N_VGetArrayPointer(y)
    Jv[0] = y[0] * v[0]
    Jv[1] = -y[0] * y[1] * v[0]
    Jv[2] = -y[1] * v[1]
    Jv[3] = y[0] * y[1] * v[1]
    return 0


def dgdu(y):
    return np.array([-1.0 + y[0], -1.0 + y[1]], dtype=np.float64)


def adjoint_rhs(t, y, l, ldot, p):
    vjp(l, ldot, t, y, p)
    ldot = N_VGetArrayPointer(ldot)
    ldot *= -1.0
    return 0


def quad_rhs(t, y, mu, qBdot, p):
    parameter_vjp(mu, qBdot, t, y, p)
    return 0


def main():
    sunctx = SUNContextView.Create()
    u = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))
    arr_u = N_VGetArrayPointer(u.get())
    arr_u[:] = 1.0

    solver = CVodeView.Create(CVodeCreate(CV_BDF, sunctx.get()))
    CVodeInit(solver.get(), lambda t, y, ydot, _: lotka_volterra(t, y, ydot, params), T0, u.get())
    CVodeSStolerances(solver.get(), RTOL, ATOL)
    ls = SUNLinearSolverView.Create(SUNLinSol_SPGMR(u.get(), 0, 3, sunctx.get()))
    CVodeSetLinearSolver(solver.get(), ls.get(), None)
    CVodeSetMaxNumSteps(solver.get(), 100000)
    CVodeAdjInit(solver.get(), STEPS, 1)  # CV_HERMITE = 1

    tout = TF
    tret = 0.0
    status, ncheck = CVodeF(solver.get(), tout, u.get(), tret, CV_NORMAL, 0)
    print(f"Forward Solution at t = {tret}: {N_VGetArrayPointer(u.get())}")

    # Adjoint terminal condition
    dg = dgdu(N_VGetArrayPointer(u.get()))
    uB = NVectorView.Create(N_VNew_Serial(NEQ, sunctx.get()))
    arr_uB = N_VGetArrayPointer(uB.get())
    arr_uB[:] = dg
    qB = NVectorView.Create(N_VNew_Serial(NP, sunctx.get()))
    arr_qB = N_VGetArrayPointer(qB.get())
    arr_qB[:] = 0.0
    print("Adjoint terminal condition:")
    print(arr_uB)
    print(arr_qB)

    # Backward problem
    which = 0
    CVodeCreateB(solver.get(), CV_BDF, which)
    CVodeInitB(
        solver.get(),
        which,
        lambda t, y, l, ldot, _: adjoint_rhs(t, y, l, ldot, params),
        TF,
        uB.get(),
    )
    CVodeSStolerancesB(solver.get(), which, RTOL, ATOL)
    lsb = SUNLinearSolverView.Create(SUNLinSol_SPGMR(uB.get(), 0, 3, sunctx.get()))
    CVodeSetLinearSolverB(solver.get(), which, lsb.get(), None)
    CVodeQuadInitB(
        solver.get(), which, lambda t, y, yB, yBdot, _: quad_rhs(t, y, yB, yBdot, params), qB.get()
    )
    CVodeSetQuadErrConB(solver.get(), which, True)
    CVodeQuadSStolerancesB(solver.get(), which, RTOL, ATOL)

    # Integrate the adjoint ODE
    CVodeB(solver.get(), T0, CV_NORMAL)
    t = 0.0
    CVodeGetB(solver.get(), which, t, uB.get())
    CVodeGetQuadB(solver.get(), which, t, qB.get())
    arr_qB = N_VGetArrayPointer(qB.get())
    arr_qB *= -1.0
    print(f"Adjoint Solution at t = {t}:")
    print(N_VGetArrayPointer(uB.get()))
    print(arr_qB)


if __name__ == "__main__":
    main()
