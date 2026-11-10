"""Compare independently integrated degenerate-FWM and gauged-coupler fields.

Run from the repository root:
    .venv/bin/python code/projection_numeric_check.py

No elliptic envelope is evaluated and no historical notebook/helper is changed.
The initial target state is the asymmetric physical example in the manuscript;
its reciprocal gauge is taken to be the identity at z=0. Matching constants and
lifted initial data are constructed at 70 digits before double-precision ODE
integration. The root with positive, largest imaginary part selects a chart
away from the real physical power orbit.
"""

import json

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import solve_ivp


def initial_data():
    with mp.workdps(70):
        fields = [mp.mpc('.9', '.2'), mp.mpc('.3', '-.1')]
        phase = mp.exp(-mp.j * mp.pi / 4)
        u = [fields[0] * phase, mp.conj(fields[1]) * phase]
        v = [-mp.conj(fields[0]) / phase, fields[1] / phase]
        a = [-mp.j * mp.mpf('.3'), -mp.j * mp.mpf('.2')]
        matrix = [[mp.j * mp.mpf('1.2'), -mp.j * mp.mpf('.25')],
                  [-mp.j * mp.mpf('.25'), mp.j * mp.mpf('.7')]]
        products = [u[j] * v[j] for j in range(2)]
        rho0 = -sum(products) / 2
        gamma = [value + rho0 for value in products]
        hamiltonian = (-sum(a[j] * products[j] for j in range(2))
                       - sum(matrix[j][k] * products[j] * products[k] / 2
                             for j in range(2) for k in range(2))
                       + mp.fprod(u) + mp.fprod(v))
        b0 = (hamiltonian + sum(a[j] * gamma[j] for j in range(2))
              + sum(matrix[j][k] * gamma[j] * gamma[k] / 2
                    for j in range(2) for k in range(2)))
        b1 = -sum(a) - sum(matrix[j][k] * gamma[k]
                            for j in range(2) for k in range(2))
        b2 = sum(value for row in matrix for value in row) / 2
        coefficients = [b0**2 - 4 * mp.fprod(gamma), 2 * b0 * b1,
                        2 * b0 * b2 + b1**2 - 4, 2 * b1 * b2, b2**2]
        roots = np.roots([complex(value) for value in coefficients[::-1]])
        seed = max(roots, key=lambda value: value.imag)
        root = mp.findroot(lambda value: mp.polyval(coefficients[::-1], value),
                           mp.mpc(complex(seed)))
        qroot = b0 + b1 * root + b2 * root**2
        d5 = sum(k * coefficients[k] * root**(k - 1) for k in range(1, 5))
        chi = 8 * qroot / d5
        factors = [mp.sqrt(2) * mp.sqrt(value - root) for value in gamma]
        sign = mp.fprod(factors) / qroot
        assert abs(sign**2 - 1) < mp.mpf('1e-60')
        sign = 1 if mp.re(sign) > 0 else -1
        theta = -b2 * gamma[0] + 2 * gamma[0] / qroot

        # Fix the ratio map's free lift scaling by bar u_3(0)=1.
        j30 = d5 / (4 * (root - rho0))
        bu = [u[j] * j30 / factors[j] for j in range(2)] + [mp.mpc(1)]
        bv = [v[j] / factors[j] for j in range(2)] + [j30]
        js = [bu[j] * bv[j] for j in range(3)]
        for j in range(2):
            assert abs(2 * js[j] - js[2] - d5 / (4 * (gamma[j] - root))) < mp.mpf('1e-60')
        energy = (-sign * chi * (bu[0] * bu[1] * bu[2]**2
                                + bv[0] * bv[1] * bv[2]**2) / 4
                  + chi * (js[0]**2 + js[1]**2) / 4 + chi * js[2]**2 / 8
                  - sum(js) / chi)
        expected_energy = -b2 * d5 / 8 - 2 / chi + 64 * (d5 + 16 * root) * root / (chi**3 * d5**2)
        assert abs(energy - expected_energy) < mp.mpf('1e-60')

        return dict(target=np.array([complex(value) for value in u + v]),
                    source=np.array([complex(value) for value in bu + bv]),
                    factors=np.array([complex(value) for value in factors]),
                    b0=complex(b0), b1=complex(b1), b2=complex(b2),
                    gamma=np.array([complex(value) for value in gamma]),
                    root=complex(root), d5=complex(d5), chi=complex(chi),
                    theta=complex(theta), sign=sign, energy=complex(energy))


def source_rhs(_z, state, chi, sign):
    u1, u2, u3, v1, v2, v3 = state
    # Direct canonical equations of the printed three-pair Hamiltonian.
    return np.array([
        -u1 / chi + chi * u1**2 * v1 / 2 - sign * chi * v2 * v3**2 / 4,
        -u2 / chi + chi * u2**2 * v2 / 2 - sign * chi * v1 * v3**2 / 4,
        -u3 / chi + chi * u3**2 * v3 / 4 - sign * chi * v1 * v2 * v3 / 2,
        v1 / chi - chi * v1**2 * u1 / 2 + sign * chi * u2 * u3**2 / 4,
        v2 / chi - chi * v2**2 * u2 / 2 + sign * chi * u1 * u3**2 / 4,
        v3 / chi - chi * v3**2 * u3 / 4 + sign * chi * u1 * u2 * u3 / 2,
    ])


def target_rhs(_z, state, b1, b2):
    u1, u2, v1, v2 = state
    return np.array([(b1 - 2 * b2 * u1 * v1) * u1 / 2 + v2,
                     (b1 - 2 * b2 * u2 * v2) * u2 / 2 + v1,
                     (-b1 + 2 * b2 * u1 * v1) * v1 / 2 - u2,
                     (-b1 + 2 * b2 * u2 * v2) * v2 / 2 - u1])


def project(z, state, data):
    u, v = state[:3], state[3:]
    phases = np.exp(data['theta'] * np.array([1, -1])[:, None] * z)
    return np.vstack((data['factors'][:, None] * u[:2] * phases / v[2],
                      data['factors'][:, None] * v[:2] / phases / u[2]))


def compare(data, rtol, atol):
    z = np.linspace(0, 8, 501)
    source = solve_ivp(source_rhs, (0, 8), data['source'],
                       args=(data['chi'], data['sign']), method='DOP853',
                       t_eval=z, rtol=rtol, atol=atol)
    target = solve_ivp(target_rhs, (0, 8), data['target'],
                       args=(data['b1'], data['b2']), method='DOP853',
                       t_eval=z, rtol=rtol, atol=atol)
    assert source.success, source.message
    assert target.success, target.message
    projected = project(z, source.y, data)
    assert np.all(np.isfinite(projected))
    error = float(np.max(np.abs(projected - target.y)))
    powers = source.y[:3] * source.y[3:]
    energy = (-data['sign'] * data['chi']
              * (source.y[0] * source.y[1] * source.y[2]**2
                 + source.y[3] * source.y[4] * source.y[5]**2) / 4
              + data['chi'] * (powers[0]**2 + powers[1]**2) / 4
              + data['chi'] * powers[2]**2 / 8 - np.sum(powers, axis=0) / data['chi'])
    differences = 2 * powers[:2] - powers[2]
    return dict(rtol=rtol, atol=atol, samples=len(z), interval=[0, 8],
                max_absolute_field_difference=error,
                scale_normalised_field_difference=error / float(np.max(np.abs(target.y))),
                max_absolute_energy_drift=float(np.max(np.abs(energy - data['energy']))),
                max_absolute_weighted_power_drift=float(np.max(np.abs(differences - differences[:, :1]))),
                sampled_min_abs_bar_u3=float(np.min(np.abs(source.y[2]))),
                sampled_min_abs_bar_v3=float(np.min(np.abs(source.y[5]))),
                source_rhs_evaluations=source.nfev, target_rhs_evaluations=target.nfev)


def main():
    data = initial_data()
    initial_error = float(np.max(np.abs(project(np.array([0.]), data['source'][:, None], data)[:, 0]
                                        - data['target'])))
    assert initial_error < 1e-14
    runs = [compare(data, 1e-11, 1e-13), compare(data, 2.5e-13, 2.5e-15)]
    assert runs[1]['max_absolute_field_difference'] < 1e-10
    assert runs[1]['max_absolute_field_difference'] < runs[0]['max_absolute_field_difference'] / 2
    assert runs[1]['max_absolute_energy_drift'] < 1e-10
    assert runs[1]['max_absolute_weighted_power_drift'] < 1e-10
    assert min(runs[1]['sampled_min_abs_bar_u3'], runs[1]['sampled_min_abs_bar_v3']) > 1e-3
    def pair(value):
        return [float(value.real), float(value.imag)]
    report = dict(description='Independent DOP853 source/target integrations; ratio projection applied only after source integration.',
                  numpy_version=np.__version__, scipy_version=scipy.__version__, mpmath_version=mp.__version__,
                  setup_precision_digits=70, b0=pair(data['b0']), b1=pair(data['b1']), b2=pair(data['b2']),
                  gamma=[pair(value) for value in data['gamma']],
                  lambda1=pair(data['root']), d5=pair(data['d5']), chi=pair(data['chi']),
                  theta=pair(data['theta']), varsigma=data['sign'],
                  source_initial_state=[pair(value) for value in data['source']],
                  target_initial_state=[pair(value) for value in data['target']],
                  initial_projection_error=initial_error, runs=runs)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
