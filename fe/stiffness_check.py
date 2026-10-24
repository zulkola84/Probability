"""Beam-theory stiffness at the loading rollers vs the digitised experiment and the paper's FE.

Shows why a plain rectangular section cannot match the paper's stress and stiffness together.
Usage: python stiffness_check.py
"""
import params as p


def delta_per_N(I, a=p.LEVER_ARM, L=p.SUPPORT_SPAN, E=p.E):
    """Deflection under a loading roller per newton of total force (two loads F/2 at distance a from the supports)."""
    def y(x, c, P):
        if x <= c:
            return P * (L - c) * x * (L**2 - (L - c) ** 2 - x**2) / (6 * E * I * L)
        return P * c * (L - x) * (L**2 - c**2 - (L - x) ** 2) / (6 * E * I * L)
    return y(a, a, 0.5) + y(a, L - a, 0.5)


if __name__ == "__main__":
    rect = p.WIDTH * p.THICKNESS**3 / 12
    print(f"Rectangle w={p.WIDTH} mm: I={rect:.1f} mm^4 -> {1 / delta_per_N(rect):.0f} N/mm (hole ignored)")
    print(f"Paper I={p.I_BENDING} mm^4        -> {1 / delta_per_N(p.I_BENDING):.0f} N/mm")
    print("Digitised experiment (Fig. 6): 363-394 N/mm;  paper FE: about 372 N/mm")
    print("The rectangle that matches the paper's FE peak stress is too stiff; the paper's own I is too flexible.")
