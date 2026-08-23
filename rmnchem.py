# RMNChem -- version MicroPython pour Casio Graph 90+E
#
# A utiliser quand on ne veut pas installer le fxSDK. C'est une version
# reduite : la base ne contient que 16 molecules et le trace est plus lent,
# casioplot n'exposant que set_pixel et draw_string.
#
# Limite importante : le MicroPython officiel de Casio n'offre aucune API
# de lecture du clavier. Les indices se saisissent donc au prompt avec
# input(), et non avec les touches de fonction comme dans la version C.
#
# Installation : copier ce fichier a la racine de la calculatrice vue en
# cle USB, puis l'ouvrir depuis le menu Python.

from casioplot import set_pixel, draw_string, show_screen, clear_screen

W, H = 384, 192            # zone de dessin de casioplot sur Graph 90+E
NOIR = (0, 0, 0)
GRIS = (130, 130, 130)
BLEU = (20, 40, 160)
ROUGE = (200, 30, 30)

# ---------------------------------------------------------------- base
# (a, b, c, d) : nom, [(delta*10, nH, multiplicite, labile), ...]
BASE = {
    (2, 6, 1, 0): [("Ethanol", [(12, 3, "t", 0), (23, 1, "s", 1),
                                (37, 2, "q", 0)]),
                   ("Dimethylether", [(33, 6, "s", 0)])],
    (3, 8, 1, 0): [("Propan-1-ol", [(9, 3, "t", 0), (16, 2, "sext", 0),
                                    (25, 1, "s", 1), (36, 2, "t", 0)]),
                   ("Propan-2-ol", [(12, 6, "d", 0), (22, 1, "s", 1),
                                    (39, 1, "sept", 0)])],
    (3, 6, 1, 0): [("Propanone", [(21, 6, "s", 0)]),
                   ("Propanal", [(11, 3, "t", 0), (24, 2, "q", 0),
                                 (98, 1, "t", 0)])],
    (2, 4, 2, 0): [("Acide ethanoique", [(21, 3, "s", 0), (114, 1, "s", 1)]),
                   ("Methanoate de methyle", [(38, 3, "s", 0),
                                              (81, 1, "s", 0)])],
    (4, 8, 2, 0): [("Ethanoate d'ethyle", [(13, 3, "t", 0), (21, 3, "s", 0),
                                           (41, 2, "q", 0)])],
    (7, 8, 0, 0): [("Toluene", [(24, 3, "s", 0), (72, 5, "m", 0)])],
    (6, 6, 0, 0): [("Benzene", [(73, 6, "s", 0)])],
    (6, 6, 1, 0): [("Phenol", [(54, 1, "s", 1), (68, 3, "m", 0),
                               (72, 2, "m", 0)])],
    (7, 6, 1, 0): [("Benzaldehyde", [(75, 3, "m", 0), (79, 2, "m", 0),
                                     (100, 1, "s", 0)])],
    (2, 5, 1, 1): [("Acetamide", [(20, 3, "s", 0), (59, 2, "s", 1)]),
                   ("N-methylformamide", [(29, 3, "d", 0), (64, 1, "s", 1),
                                          (81, 1, "s", 0)])],
    (3, 7, 2, 1): [("Alanine", [(15, 3, "d", 0), (38, 1, "q", 0),
                                (60, 3, "s", 1)])],
}

RAIES = {"s": 1, "d": 2, "t": 3, "q": 4, "quint": 5, "sext": 6,
         "sept": 7, "m": 3}


# ------------------------------------------------------------- dessin
def ligne(x1, y1, x2, y2, c):
    """Bresenham : casioplot ne fournit pas de primitive de ligne."""
    x1, y1, x2, y2 = int(x1), int(y1), int(x2), int(y2)
    dx = abs(x2 - x1)
    dy = -abs(y2 - y1)
    sx = 1 if x1 < x2 else -1
    sy = 1 if y1 < y2 else -1
    err = dx + dy
    while True:
        if 0 <= x1 < W and 0 <= y1 < H:
            set_pixel(x1, y1, c)
        if x1 == x2 and y1 == y2:
            return
        e2 = 2 * err
        if e2 >= dy:
            err += dy
            x1 += sx
        if e2 <= dx:
            err += dx
            y1 += sy


def formule(a, b, c, d):
    s = ""
    for sym, n in (("C", a), ("H", b), ("O", c), ("N", d)):
        if n:
            s += sym + (str(n) if n > 1 else "")
    return s or "-"


def degre_insaturation(a, b, c, d):
    """Renvoie (2*DI, message d'erreur eventuel)."""
    if a == b == c == d == 0:
        return None, "saisir au moins un atome"
    if b > 2 * a + 2 + d:
        return None, "trop de H : b depasse 2a+2+d"
    di2 = 2 * a + 2 + d - b
    if di2 < 0:
        return None, "DI negatif"
    if di2 % 2:
        return None, "DI demi-entier : ion ou radical"
    return di2, None


def fragments(a, b, c, d, di):
    out = []
    if di == 0:
        out.append("chaine saturee acyclique")
    if di >= 4 and a >= 6:
        out.append("noyau benzenique probable")
    if 1 <= di <= 3:
        out.append("C=C, C=O ou cycle")
    if c == 1:
        out.append("C=O aldehyde/cetone" if di >= 1 else "alcool ou ether")
    if c >= 2:
        out.append("acide ou ester" if di >= 1 else "diol ou diether")
    if d >= 1 and c == 0:
        out.append("nitrile (C#N vaut 2)" if di >= 2 else "amine")
    if d >= 1 and c >= 1 and di >= 1:
        out.append("amide ou derive nitre")
    return out


PPM_MAX = 120
SX0, SX1 = 24, W - 8
SBASE = H - 24


def x_ppm(s10):
    s10 = max(0, min(PPM_MAX, s10))
    return SX1 - (SX1 - SX0) * s10 // PPM_MAX


def trace_spectre(pics):
    ligne(SX0 - 4, SBASE, SX1 + 2, SBASE, NOIR)
    for p in range(0, 13, 2):
        x = x_ppm(p * 10)
        ligne(x, SBASE, x, SBASE + 4, GRIS)
        draw_string(x - 4, SBASE + 6, str(p), GRIS, "small")
    draw_string(0, SBASE + 6, "ppm", GRIS, "small")

    hmax = max(p[1] for p in pics)
    for delta, nH, mult, labile in pics:
        x = x_ppm(delta)
        h = 14 + nH * 46 // hmax
        n = RAIES[mult]
        col = BLEU if not labile else (60, 120, 200)
        for k in range(n):
            xx = x - (n - 1) + 2 * k
            hh = h - h * abs(k - (n - 1) // 2) // (n + 2)
            ligne(xx, SBASE - 1, xx, SBASE - hh, col)
        draw_string(x - 8, SBASE - h - 16, "%dH" % nH, col, "small")
        draw_string(x - 8, SBASE - h - 30, mult, GRIS, "small")


def affiche(nom, a, b, c, d, di, pics):
    clear_screen()
    draw_string(2, 2, formule(a, b, c, d), NOIR, "medium")
    draw_string(110, 2, nom, BLEU, "medium")
    draw_string(300, 2, "DI=%d" % di, ROUGE, "medium")
    y = 26
    for delta, nH, mult, labile in sorted(pics):
        txt = "%d,%d ppm  %dH  %s%s" % (delta // 10, delta % 10, nH, mult,
                                        "  (labile)" if labile else "")
        draw_string(4, y, txt, NOIR, "small")
        y += 14
    trace_spectre(pics)
    show_screen()


def secours(a, b, c, d, di, err):
    clear_screen()
    draw_string(2, 2, formule(a, b, c, d), NOIR, "medium")
    if err:
        draw_string(4, 26, "Formule impossible", ROUGE, "medium")
        draw_string(4, 48, err, NOIR, "small")
        show_screen()
        return
    draw_string(200, 2, "DI = %d" % di, ROUGE, "medium")
    draw_string(4, 26, "Non tabulee. Fragments :", BLEU, "small")
    y = 44
    for f in fragments(a, b, c, d, di):
        draw_string(10, y, "- " + f, NOIR, "small")
        y += 14
    draw_string(10, y + 6, "H a placer : %d" % b, NOIR, "small")
    trace_spectre([(15, 1, "s", 0)])
    show_screen()


def demande(nom, defaut):
    try:
        s = input("%s [%d] : " % (nom, defaut)).strip()
        return int(s) if s else defaut
    except (ValueError, EOFError):
        return defaut


def main():
    print("RMNChem -- formule CaHbOcNd")
    print("Entree vide = valeur par defaut, 0 = element absent")
    a = demande("a (C)", 2)
    b = demande("b (H)", 6)
    c = demande("c (O)", 1)
    d = demande("d (N)", 0)

    di2, err = degre_insaturation(a, b, c, d)
    di = di2 // 2 if di2 is not None else 0

    if err:
        secours(a, b, c, d, di, err)
        return

    trouves = BASE.get((a, b, c, d))
    if not trouves:
        secours(a, b, c, d, di, None)
        return

    for i, (nom, pics) in enumerate(trouves):
        affiche(nom, a, b, c, d, di, pics)
        if i < len(trouves) - 1:
            input("isomere %d/%d -- EXE pour le suivant "
                  % (i + 1, len(trouves)))


main()
