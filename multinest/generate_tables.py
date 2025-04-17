import numpy as np
from typing import Callable, TYPE_CHECKING
import numpy.typing
from scipy.fft import dct
from scipy.integrate import quad as nintegrate
from scipy.optimize import fsolve


nda = numpy.typing.NDArray[np.float64]

if TYPE_CHECKING:

    def exp(x: float) -> float:
        return np.exp(x)

    def sqrt(x: float) -> float:
        return np.sqrt(x)

    def cos(x: nda) -> nda:
        return np.cos(x)

else:
    from numpy import cos, exp, sqrt


def nodes(n: int) -> nda:
    x = cos(np.arange(0.5, n + 0.5 + 1) * np.pi / (n + 1))
    return x


def coeff(n: int, func: Callable[[float], float]) -> nda:
    vals = np.array([func(x) for x in nodes(n)])
    coeffs = 2 * dct(vals, norm="forward")
    assert isinstance(coeffs, np.ndarray)
    coeffs[0] /= 2
    return coeffs


def build_lvc_int(n: int = 30) -> nda:
    def integrand(x: float) -> float:
        return 1 / (1 + exp((1 - 2 * x) / (x - x**2)))

    return coeff(
        n,
        lambda z: nintegrate(integrand, 0, (1 + z) / 2, epsabs=1e-10, epsrel=1e-10)[0],
    )


def build_H0(luminosity_relation: Callable[[float], float], n: int = 100) -> nda:
    conv = 299792.458

    def Dl(z: float) -> float:
        return (
            conv
            * (1 + z)
            * nintegrate(luminosity_relation, 0, z, epsabs=1e-10, epsrel=1e-10)[0]
        )

    def func(y: float) -> float:
        return fsolve(lambda z: y - Dl(z) / conv / 10 + 1, 1)[0]

    return coeff(n, func)


def build_H0_Planck(n: int = 100, Om: float = 0.3111) -> nda:
    return build_H0(lambda zz: 1 / sqrt(Om * (1 + zz) ** 3 + (1 - Om)), n)


def to_fortran(dat: nda, n: int = 3) -> str:
    body = ", &\n".join(
        "      "
        + ", ".join(f"{'+' if i>0 else '-'}{abs(i):.15e}_prec" for i in dat[j : j + n])
        for j in range(0, len(dat), n)
    )
    return f"  integer, parameter :: n = {len(dat)}\n  real(kind=prec), parameter :: tab(n) =  (/&\n{body} /)"


if __name__ == "__main__":
    import sys

    if sys.argv[1] == "lvc":
        dat = build_lvc_int()
    elif sys.argv[1] == "Planck":
        dat = build_H0_Planck()

    print(to_fortran(dat))
