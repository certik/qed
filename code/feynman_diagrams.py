"""Generate SVG Feynman diagrams for the g-2 notes (figures/*.svg).

Geometry: the electron comes in from the lower left, goes up to the apex
(where the external photon gamma(q) attaches, drawn as a vertical wavy
line), and leaves to the lower right.  Vertex positions on the legs are
given as fractions: s in [0,1] along the incoming leg (0 = start,
1 = apex), r in [0,1] along the outgoing leg (0 = apex, 1 = end).
"""
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

APEX = np.array([3.0, 3.0])
IN0 = np.array([0.0, 0.0])
OUT1 = np.array([6.0, 0.0])


def leg_point(leg, f):
    if leg == "in":
        return IN0 + f * (APEX - IN0)
    return APEX + f * (OUT1 - APEX)


def wavy_path(p0, p1, amp=0.10, halfwaves=8, n=400):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    L = np.hypot(*d)
    t = np.linspace(0, 1, n)
    perp = np.array([-d[1], d[0]]) / L
    ripple = amp * np.sin(np.pi * halfwaves * t) * np.sin(np.pi * t) ** 0.1
    pts = p0[None, :] + t[:, None] * d[None, :] + ripple[:, None] * perp[None, :]
    return pts


def wavy_arc(pa, pb, amp=0.08, halfwaves=8, bulge=1.0, n=400):
    """Wavy line along a circular-ish arc bulging to the left of pa->pb."""
    pa, pb = np.asarray(pa, float), np.asarray(pb, float)
    d = pb - pa
    L = np.hypot(*d)
    perp = np.array([-d[1], d[0]]) / L
    t = np.linspace(0, 1, n)
    base = pa[None, :] + t[:, None] * d[None, :]
    arc = bulge * 0.5 * L * np.sin(np.pi * t)
    ripple = amp * np.sin(np.pi * halfwaves * t) * np.sin(np.pi * t) ** 0.1
    pts = base + (arc + ripple)[:, None] * perp[None, :]
    return pts


def draw_electron_line(ax):
    for a, b in [(IN0, APEX), (APEX, OUT1)]:
        ax.plot([a[0], b[0]], [a[1], b[1]], "k-", lw=1.6, zorder=2)
        mid = 0.5 * (a + b)
        d = (b - a) / np.hypot(*(b - a))
        ax.annotate(
            "",
            xy=mid + 0.12 * d,
            xytext=mid - 0.12 * d,
            arrowprops=dict(arrowstyle="-|>", color="k", lw=1.6),
            zorder=3,
        )


def draw_external_photon(ax):
    pts = wavy_path(APEX, APEX + np.array([0, 1.5]), amp=0.09, halfwaves=6)
    ax.plot(pts[:, 0], pts[:, 1], "k-", lw=1.2)
    ax.text(APEX[0] + 0.22, APEX[1] + 1.25, r"$\gamma(q)$", fontsize=13)


def draw_photon(ax, a, b, **kw):
    pts = wavy_path(a, b, **kw)
    ax.plot(pts[:, 0], pts[:, 1], "k-", lw=1.2, zorder=1)


def draw_photon_with_bubble(ax, a, b, rad=0.42):
    a, b = np.asarray(a, float), np.asarray(b, float)
    mid = 0.5 * (a + b)
    d = (b - a) / np.hypot(*(b - a))
    c1, c2 = mid - rad * d, mid + rad * d
    draw_photon(ax, a, c1, halfwaves=5)
    draw_photon(ax, c2, b, halfwaves=5)
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(mid[0] + rad * np.cos(th), mid[1] + rad * np.sin(th), "k-", lw=1.6)
    for ang in (np.pi / 2, 3 * np.pi / 2):
        p = mid + rad * np.array([np.cos(ang), np.sin(ang)])
        d2 = np.array([-np.sin(ang), np.cos(ang)])
        ax.annotate(
            "",
            xy=p + 0.1 * d2,
            xytext=p - 0.1 * d2,
            arrowprops=dict(arrowstyle="-|>", color="k", lw=1.4),
        )


def draw_se_arc(ax, leg, f1, f2):
    """Photon arc (self-energy) between two points of the same leg."""
    a, b = leg_point(leg, f1), leg_point(leg, f2)
    if leg == "in":
        a, b = b, a  # bulge away from the diagram interior
    pts = wavy_arc(a, b, bulge=0.8)
    ax.plot(pts[:, 0], pts[:, 1], "k-", lw=1.2, zorder=1)


def draw_cross(ax, p, size=0.14):
    ax.plot(
        [p[0] - size, p[0] + size], [p[1] - size, p[1] + size], "k-", lw=2.2, zorder=4
    )
    ax.plot(
        [p[0] - size, p[0] + size], [p[1] + size, p[1] - size], "k-", lw=2.2, zorder=4
    )


def new_fig():
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    ax.set_xlim(-0.7, 6.7)
    ax.set_ylim(-0.7, 4.9)
    ax.set_aspect("equal")
    ax.axis("off")
    draw_electron_line(ax)
    draw_external_photon(ax)
    ax.text(-0.55, -0.05, r"$p$", fontsize=13)
    ax.text(6.25, -0.05, r"$p'$", fontsize=13)
    return fig, ax


def save(fig, name):
    outdir = os.path.join(os.path.dirname(__file__), "..", "figures")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.normpath(path))


def fig_lo():
    fig, ax = new_fig()
    draw_photon(ax, leg_point("in", 0.4), leg_point("out", 0.6))
    save(fig, "g2-lo.svg")


def fig_I():  # crossed ladder
    fig, ax = new_fig()
    draw_photon(ax, leg_point("in", 0.3), leg_point("out", 0.35), halfwaves=10)
    draw_photon(ax, leg_point("in", 0.65), leg_point("out", 0.7), halfwaves=10)
    save(fig, "g2-nlo-I.svg")


def fig_IIa():  # ladder (nested)
    fig, ax = new_fig()
    draw_photon(ax, leg_point("in", 0.3), leg_point("out", 0.7), halfwaves=12)
    draw_photon(ax, leg_point("in", 0.65), leg_point("out", 0.35), halfwaves=8)
    save(fig, "g2-nlo-IIa.svg")


def fig_IIb():  # self-energy on external leg
    fig, ax = new_fig()
    draw_photon(ax, leg_point("in", 0.55), leg_point("out", 0.5))
    draw_se_arc(ax, "in", 0.08, 0.35)
    save(fig, "g2-nlo-IIb.svg")


def fig_IIc():  # corner: vertex part at an internal vertex
    fig, ax = new_fig()
    draw_photon(ax, leg_point("in", 0.4), leg_point("out", 0.55), halfwaves=10)
    draw_se_arc(ax, "out", 0.3, 0.8)
    save(fig, "g2-nlo-IIc.svg")


def fig_IId():  # self-energy insertion on internal line
    fig, ax = new_fig()
    draw_photon(ax, leg_point("in", 0.4), leg_point("out", 0.8), halfwaves=12)
    draw_se_arc(ax, "out", 0.3, 0.6)
    save(fig, "g2-nlo-IId.svg")


def fig_IIe():  # vacuum polarization insertion
    fig, ax = new_fig()
    draw_photon_with_bubble(ax, leg_point("in", 0.35), leg_point("out", 0.65))
    save(fig, "g2-nlo-IIe.svg")


def fig_deltam():  # mass counterterm insertion
    fig, ax = new_fig()
    draw_photon(ax, leg_point("in", 0.4), leg_point("out", 0.6))
    draw_cross(ax, leg_point("in", 0.7))
    ax.text(
        leg_point("in", 0.7)[0] - 0.75,
        leg_point("in", 0.7)[1] + 0.05,
        r"$\delta m$",
        fontsize=13,
    )
    save(fig, "g2-nlo-deltam.svg")


if __name__ == "__main__":
    fig_lo()
    fig_I()
    fig_IIa()
    fig_IIb()
    fig_IIc()
    fig_IId()
    fig_IIe()
    fig_deltam()
