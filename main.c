/*
 * RMNChem -- analyse de formule brute CaHbOcNd
 * Casio Graph 90+E / fx-CG 50 -- gint / fxSDK
 *
 * Saisie des indices a, b, c, d -> degre d'insaturation, identification des
 * isomeres tabules, trace de la structure 2D et du spectre RMN 1H simule.
 *
 * Compilation : fxsdk build-cg
 */

#include <gint/display.h>
#include <gint/keyboard.h>

#include "chem_db.h"

/* ================================================================== */
/* Geometrie de l'ecran                                                */
/*                                                                     */
/* gint donne acces a la dalle entiere : 396 x 224 (DWIDTH x DHEIGHT). */
/* La zone 384 x 216 souvent citee est celle que l'OS Casio laisse aux */
/* programmes Basic ; sous gint la bordure est utilisable. Tout est    */
/* exprime a partir de DWIDTH / DHEIGHT, donc le code suit la dalle    */
/* reelle quelle que soit la machine de la famille fx-CG.              */
/* ================================================================== */
#define W          DWIDTH
#define H          DHEIGHT

#define CHW        8              /* largeur d'un caractere            */
#define CHH        9              /* hauteur d'un caractere            */

#define HDR_H      26             /* bandeau superieur                 */
#define FKEY_Y     (H - 18)       /* bandeau de touches de fonction    */

#define STRUCT_X0  2
#define STRUCT_X1  172
#define TABLE_X0   176
#define TABLE_X1   (W - 4)
#define MID_Y0     (HDR_H + 4)
#define MID_Y1     116

#define SPEC_X0    30
#define SPEC_X1    (W - 6)
#define SPEC_TOP   122
#define SPEC_BASE  (FKEY_Y - 22)  /* ligne de base du spectre          */

#define FB_X       208            /* colonne droite du mode secours    */
#define PPM_MAX    120            /* 12,0 ppm : borne gauche de l'axe  */

/* Couleurs */
#define C_GREY     C_RGB(20, 20, 20)
#define C_DARKB    C_RGB(0, 0, 18)
#define C_ORANGE   C_RGB(31, 18, 0)

/* ================================================================== */
/* Etat global                                                         */
/* ================================================================== */
static uint8_t in_a = 2, in_b = 6, in_c = 1, in_d = 0;
static int8_t  cur_field = 0;      /* 0=a 1=b 2=c 3=d                  */
static int8_t  typing = 0;         /* saisie numerique en cours        */

static int16_t hits[16];           /* index des isomeres trouves       */
static int8_t  nhits, cur_hit;
static int8_t  show_integral = 1;

/* ================================================================== */
/* Utilitaires                                                         */
/*                                                                     */
/* Pas de flottant ici : dprint de gint ne formate pas les reels et le */
/* SH4 de la Graph 90+E n'a pas d'unite double precision. Les          */
/* deplacements chimiques circulent en dixiemes de ppm.                */
/* ================================================================== */

static int slen(const char *s)
{
	int n = 0;
	while(s[n]) n++;
	return n;
}

static int u2s(char *p, unsigned v)
{
	char tmp[6];
	int n = 0, i;
	if(!v) { *p = '0'; return 1; }
	while(v) { tmp[n++] = (char)('0' + v % 10); v /= 10; }
	for(i = 0; i < n; i++) p[i] = tmp[n - 1 - i];
	return n;
}

/* 72 -> "7,2" (virgule decimale) */
static int fmt_shift(char *p, unsigned v10)
{
	int n = u2s(p, v10 / 10);
	p[n++] = ',';
	p[n++] = (char)('0' + (v10 % 10));
	return n;
}

/* Formule brute avec omission des indices 0 et 1 : C6H4O2 */
static void fmt_formula(char *buf, int a, int b, int c, int d)
{
	int n = 0;
	if(a) { buf[n++] = 'C'; if(a > 1) n += u2s(buf + n, a); }
	if(b) { buf[n++] = 'H'; if(b > 1) n += u2s(buf + n, b); }
	if(c) { buf[n++] = 'O'; if(c > 1) n += u2s(buf + n, c); }
	if(d) { buf[n++] = 'N'; if(d > 1) n += u2s(buf + n, d); }
	if(!n) buf[n++] = '-';
	buf[n] = 0;
}

static void atom_label(char *buf, const atom_t *at)
{
	int n = 0;
	buf[n++] = (at->el == EL_C) ? 'C' : (at->el == EL_O) ? 'O' : 'N';
	if(at->nH) {
		buf[n++] = 'H';
		if(at->nH > 1) buf[n++] = (char)('0' + at->nH);
	}
	buf[n] = 0;
}

/* Racine carree entiere : sert a raccourcir les liaisons sous les labels */
static int isqrt32(int v)
{
	int r = 0, b = 1 << 14;
	if(v <= 0) return 0;
	while(b > v) b >>= 2;
	while(b) {
		if(v >= r + b) { v -= r + b; r = (r >> 1) + b; }
		else r >>= 1;
		b >>= 2;
	}
	return r;
}

/* ================================================================== */
/* Moteur : degre d'insaturation                                       */
/*                                                                     */
/*   DI = a + 1 + (d - b) / 2                                          */
/*                                                                     */
/* On calcule 2*DI pour rester en entiers : une valeur impaire signale */
/* un DI demi-entier, donc une formule qui ne peut pas correspondre a  */
/* une molecule neutre a couches completes.                            */
/* ================================================================== */
enum { DI_OK = 0, DI_HALF, DI_NEG, DI_TOOMANY_H, DI_EMPTY };

static int di_x2;

static int compute_di(int a, int b, int c, int d)
{
	(void)c;
	if(!a && !b && !c && !d) return DI_EMPTY;
	di_x2 = 2 * a + 2 + d - b;
	if(b > 2 * a + 2 + d) return DI_TOOMANY_H;
	if(di_x2 < 0) return DI_NEG;
	if(di_x2 & 1) return DI_HALF;
	return DI_OK;
}

static void search_db(int a, int b, int c, int d)
{
	int i;
	nhits = 0;
	cur_hit = 0;
	for(i = 0; i < DB_SIZE && nhits < 16; i++)
		if(DB[i].a == a && DB[i].b == b && DB[i].c == c && DB[i].d == d)
			hits[(int)nhits++] = (int16_t)i;
}

/* ================================================================== */
/* Trace de la structure                                               */
/* ================================================================== */

static void circle_outline(int cx, int cy, int r, int color)
{
	int x = r, y = 0, err = 1 - r;
	while(x >= y) {
		dpixel(cx + x, cy + y, color); dpixel(cx + y, cy + x, color);
		dpixel(cx - y, cy + x, color); dpixel(cx - x, cy + y, color);
		dpixel(cx - x, cy - y, color); dpixel(cx - y, cy - x, color);
		dpixel(cx + y, cy - x, color); dpixel(cx + x, cy - y, color);
		y++;
		if(err < 0) err += 2 * y + 1;
		else { x--; err += 2 * (y - x) + 1; }
	}
}

/*
 * Liaison, eventuellement multiple, raccourcie a ses extremites pour
 * laisser la place aux etiquettes d'atomes. Le decalage des traits
 * supplementaires approxime la normale par l'axe dominant : pas de
 * calcul en plus, et le rendu reste correct sur les geometries de la
 * base, qui n'utilisent que des directions franches.
 */
static void draw_bond(int x1, int y1, int x2, int y2, int order,
                      int cut1, int cut2, int color)
{
	int dx = x2 - x1, dy = y2 - y1;
	int len = isqrt32(dx * dx + dy * dy);
	int ox, oy, k, adx = dx < 0 ? -dx : dx, ady = dy < 0 ? -dy : dy;

	if(len < 4 || cut1 + cut2 >= len - 2) return;

	x1 += dx * cut1 / len; y1 += dy * cut1 / len;
	x2 -= dx * cut2 / len; y2 -= dy * cut2 / len;

	if(adx >= ady) { ox = 0; oy = 3; } else { ox = 3; oy = 0; }

	dline(x1, y1, x2, y2, color);
	if(order == 2)
		dline(x1 + ox, y1 + oy, x2 + ox, y2 + oy, color);
	else if(order == 3)
		for(k = -1; k <= 1; k += 2)
			dline(x1 + k * ox, y1 + k * oy,
			      x2 + k * ox, y2 + k * oy, color);
}

/*
 * Rendu d'une molecule. La boite englobante est mise a l'echelle pour
 * remplir le cadre, avec un facteur exprime en huitiemes et borne entre
 * 1,0 et 2,5 : les coordonnees de la base restent ainsi de simples
 * entiers relatifs, lisibles et compacts.
 */
static void draw_molecule(const molecule_t *m, int rx0, int ry0,
                          int rx1, int ry1)
{
	int i, minx = 127, maxx = -128, miny = 127, maxy = -128;
	int sx, sy, s8, ox, oy, arx = 0, ary = 0, nar = 0;
	int bw, bh, pw = rx1 - rx0, ph = ry1 - ry0;
	char lbl[6];

	for(i = 0; i < m->natoms; i++) {
		const atom_t *at = &m->atoms[i];
		if(at->x < minx) minx = at->x;
		if(at->x > maxx) maxx = at->x;
		if(at->y < miny) miny = at->y;
		if(at->y > maxy) maxy = at->y;
	}
	bw = maxx - minx;
	bh = maxy - miny;

	/* on reserve de la place pour les etiquettes de bord */
	sx = bw > 0 ? (8 * (pw - 40)) / bw : 20;
	sy = bh > 0 ? (8 * (ph - 20)) / bh : 20;
	s8 = sx < sy ? sx : sy;
	if(s8 > 20) s8 = 20;
	if(s8 < 8)  s8 = 8;

	ox = (rx0 + rx1) / 2 - ((minx + maxx) / 2) * s8 / 8;
	oy = (ry0 + ry1) / 2 - ((miny + maxy) / 2) * s8 / 8;

	for(i = 0; i < m->nbonds; i++) {
		const bond_t *b = &m->bonds[i];
		const atom_t *A, *B;
		int ax, ay, bx, by, c1 = 0, c2 = 0, col, horiz;
		if(b->i == b->j) continue;                 /* liaison factice */
		if(b->i >= m->natoms || b->j >= m->natoms) continue;
		A = &m->atoms[b->i];
		B = &m->atoms[b->j];
		ax = ox + A->x * s8 / 8; ay = oy + A->y * s8 / 8;
		bx = ox + B->x * s8 / 8; by = oy + B->y * s8 / 8;
		horiz = ((bx - ax) > 10 || (ax - bx) > 10);

		if(A->show) {
			atom_label(lbl, A);
			c1 = horiz ? slen(lbl) * CHW / 2 + 3 : CHH / 2 + 3;
		}
		if(B->show) {
			atom_label(lbl, B);
			c2 = horiz ? slen(lbl) * CHW / 2 + 3 : CHH / 2 + 3;
		}
		col = (A->hl || B->hl) ? C_RED : C_BLACK;
		draw_bond(ax, ay, bx, by, b->order, c1, c2, col);
		if(b->arom) { arx += A->x; ary += A->y; nar++; }
	}

	if(nar >= 6)
		circle_outline(ox + (arx / nar) * s8 / 8,
		               oy + (ary / nar) * s8 / 8,
		               13 * s8 / 8, C_BLACK);

	for(i = 0; i < m->natoms; i++) {
		const atom_t *at = &m->atoms[i];
		if(!at->show) continue;
		atom_label(lbl, at);
		dtext_opt(ox + at->x * s8 / 8, oy + at->y * s8 / 8,
		          at->hl ? C_RED : C_BLACK, C_NONE,
		          DTEXT_CENTER, DTEXT_MIDDLE, lbl, -1);
	}
}

/* ================================================================== */
/* Spectre RMN 1H                                                      */
/* ================================================================== */

static int ppm_to_x(int shift10)
{
	if(shift10 > PPM_MAX) shift10 = PPM_MAX;
	if(shift10 < 0) shift10 = 0;
	return SPEC_X1 - ((SPEC_X1 - SPEC_X0) * shift10) / PPM_MAX;
}

/* Tri des signaux par deplacement croissant : la base reste ecrite dans
 * l'ordre de lecture de la molecule, l'affichage se charge du reste. */
static int sort_peaks(const nmr_peak_t *pk, int n, int *order)
{
	int i, j;
	if(n > 12) n = 12;
	for(i = 0; i < n; i++) order[i] = i;
	for(i = 1; i < n; i++) {
		int t = order[i];
		for(j = i; j > 0 && pk[order[j - 1]].shift > pk[t].shift; j--)
			order[j] = order[j - 1];
		order[j] = t;
	}
	return n;
}

static void draw_axis(int with_tms)
{
	char buf[4];
	int ppm;
	dline(SPEC_X0 - 4, SPEC_BASE, SPEC_X1 + 2, SPEC_BASE, C_BLACK);
	for(ppm = 0; ppm <= 12; ppm++) {
		int x = ppm_to_x(ppm * 10);
		dline(x, SPEC_BASE, x, SPEC_BASE + ((ppm % 2) ? 2 : 4), C_GREY);
		if(ppm % 2 == 0) {
			int n = u2s(buf, ppm);
			buf[n] = 0;
			dtext_opt(x, SPEC_BASE + 11, C_GREY, C_NONE,
			          DTEXT_CENTER, DTEXT_MIDDLE, buf, -1);
		}
	}
	dtext(2, SPEC_BASE + 6, C_GREY, "ppm");
	if(with_tms) {
		dline(SPEC_X1, SPEC_BASE, SPEC_X1, SPEC_BASE - 8, C_GREY);
		dtext_opt(SPEC_X1, SPEC_BASE - 14, C_GREY, C_NONE,
		          DTEXT_RIGHT, DTEXT_MIDDLE, "TMS", -1);
	}
}

static void draw_spectrum(const nmr_peak_t *pk, int n)
{
	int order[12], i, k, maxH = 1, total = 0, acc = 0;
	int prev_x = SPEC_X1, prev_y = SPEC_BASE - 3;
	int hmax = SPEC_BASE - SPEC_TOP - 20;
	char buf[6];

	n = sort_peaks(pk, n, order);
	draw_axis(1);

	for(i = 0; i < n; i++) {
		if(pk[i].nH > maxH) maxH = pk[i].nH;
		total += pk[i].nH;
	}

	/* Courbe d'integration : escalier parcouru de 0 ppm vers les champs
	 * faibles, donc de la droite vers la gauche de l'ecran. Tracee avant
	 * les pics pour rester en arriere-plan. */
	if(show_integral && total) {
		for(i = 0; i < n; i++) {
			const nmr_peak_t *p = &pk[order[i]];
			int x = ppm_to_x(p->shift), y;
			acc += p->nH;
			y = SPEC_BASE - 3 - (acc * hmax) / total;
			dline(prev_x, prev_y, x + 4, prev_y, C_ORANGE);
			dline(x + 4, prev_y, x - 4, y, C_ORANGE);
			prev_x = x - 4;
			prev_y = y;
		}
		dline(prev_x, prev_y, SPEC_X0 - 4, prev_y, C_ORANGE);
	}

	for(i = 0; i < n; i++) {
		const nmr_peak_t *p = &pk[order[i]];
		int x = ppm_to_x(p->shift);
		int nl = MULT_LINES[p->mult];
		int h = 12 + (p->nH * (hmax - 12)) / maxH;
		int col = p->exch ? C_BLUE : C_DARKB;
		int x0 = x - (nl - 1);

		for(k = 0; k < nl; k++) {
			int xx = x0 + 2 * k, hh = h;
			if(nl > 1) {
				int dc = k - (nl - 1) / 2;
				if(dc < 0) dc = -dc;
				hh = h - (h * dc) / (nl + 2);
			}
			if(p->mult == M_M) hh = h - 3 * k;
			dline(xx, SPEC_BASE - 1, xx, SPEC_BASE - hh, col);
		}
		{
			int nn = u2s(buf, p->nH);
			buf[nn++] = 'H';
			buf[nn] = 0;
			dtext_opt(x, SPEC_BASE - h - 6, col, C_NONE,
			          DTEXT_CENTER, DTEXT_MIDDLE, buf, -1);
			dtext_opt(x, SPEC_BASE - h - 15, C_GREY, C_NONE,
			          DTEXT_CENTER, DTEXT_MIDDLE,
			          MULT_TXT[p->mult], -1);
		}
	}
}

/* Axe seul avec les zones de reference, pour le moteur de secours */
static void draw_zone_axis(void)
{
	static const uint8_t z0[] = {  5, 20, 33, 45, 65,  95 };
	static const uint8_t z1[] = { 18, 30, 43, 60, 85, 105 };
	static const char * const zn[] = {
		"CH alcane", "CH-C=O", "CH-O", "OH / NH", "H arom.", "CHO"
	};
	int i;

	draw_axis(0);
	for(i = 0; i < 6; i++) {
		int xa = ppm_to_x(z1[i]), xb = ppm_to_x(z0[i]);
		int y = SPEC_BASE - 12 - 14 * (i % 3);
		int wlab = slen(zn[i]) * CHW, xc = (xa + xb) / 2;
		drect(xa, y, xb, y + 4, (i & 1) ? C_RGB(14, 22, 31)
		                                : C_RGB(20, 29, 20));
		if(xc - wlab / 2 < 2) xc = wlab / 2 + 2;
		if(xc + wlab / 2 > W - 2) xc = W - 2 - wlab / 2;
		dtext_opt(xc, y - 5, C_GREY, C_NONE,
		          DTEXT_CENTER, DTEXT_MIDDLE, zn[i], -1);
	}
}

/* ================================================================== */
/* Tableau des signaux                                                 */
/*                                                                     */
/* 27 caracteres disponibles : delta(5) nH(4) multiplicite(5) puis     */
/* l'attribution. Un asterisque marque les protons labiles.            */
/* ================================================================== */
static void draw_table(const nmr_peak_t *pk, int n)
{
	int order[12], i, y = MID_Y0 + 13, labile = 0;
	char buf[32];

	dtext(TABLE_X0, MID_Y0, C_BLUE, "ppm  nH   mult attribution");
	dline(TABLE_X0, MID_Y0 + 10, TABLE_X1, MID_Y0 + 10, C_GREY);

	n = sort_peaks(pk, n, order);
	for(i = 0; i < n; i++) {
		const nmr_peak_t *p = &pk[order[i]];
		const char *s;
		int k = fmt_shift(buf, p->shift);
		while(k < 5) buf[k++] = ' ';
		k += u2s(buf + k, p->nH);
		buf[k++] = 'H';
		if(p->exch) { buf[k++] = '*'; labile = 1; }
		while(k < 9) buf[k++] = ' ';
		for(s = MULT_TXT[p->mult]; *s; s++) buf[k++] = *s;
		while(k < 14) buf[k++] = ' ';
		for(s = ATTR[p->attr]; *s && k < 27; s++) buf[k++] = *s;
		buf[k] = 0;
		dtext(TABLE_X0, y, p->exch ? C_BLUE : C_BLACK, buf);
		y += 11;
		if(y > MID_Y1 - 20) break;
	}
	if(labile) {
		dline(TABLE_X0, MID_Y1 - 16, TABLE_X1, MID_Y1 - 16, C_GREY);
		dtext(TABLE_X0, MID_Y1 - 13, C_BLUE,
		      "* labile : efface par D2O");
	}
}

/* ================================================================== */
/* Bandeaux                                                            */
/* ================================================================== */
static void draw_fkeys(const char *f1, const char *f2, const char *f3,
                       const char *f4, const char *f5, const char *f6)
{
	const char *lab[6] = { f1, f2, f3, f4, f5, f6 };
	int i;
	drect(0, FKEY_Y, W - 1, H - 1, C_RGB(28, 28, 28));
	for(i = 0; i < 6; i++) {
		int x0 = (W * i) / 6, x1 = (W * (i + 1)) / 6 - 2;
		if(!lab[i]) continue;
		drect(x0 + 1, FKEY_Y + 2, x1, H - 3, C_WHITE);
		dtext_opt((x0 + x1) / 2, FKEY_Y + 9, C_BLACK, C_NONE,
		          DTEXT_CENTER, DTEXT_MIDDLE, lab[i], -1);
	}
}

static void draw_header(int a, int b, int c, int d, const char *name,
                        const char *sub)
{
	char buf[24];
	fmt_formula(buf, a, b, c, d);
	drect(0, 0, W - 1, HDR_H - 3, C_RGB(0, 0, 12));
	dtext(4, 2, C_WHITE, buf);
	if(name) dtext(96, 2, C_RGB(31, 31, 12), name);
	if(sub)  dtext(4, 13, C_RGB(24, 24, 24), sub);
	dline(0, HDR_H - 1, W - 1, HDR_H - 1, C_BLACK);
}

/* ================================================================== */
/* Ecran de saisie                                                     */
/* ================================================================== */
static void screen_input(void)
{
	static const char *lname[4] = { "C  (a)", "H  (b)", "O  (c)", "N  (d)" };
	uint8_t *val[4] = { &in_a, &in_b, &in_c, &in_d };
	char buf[24];
	int i, st;

	dclear(C_WHITE);
	drect(0, 0, W - 1, 26, C_RGB(0, 0, 14));
	dtext_opt(W / 2, 7, C_WHITE, C_NONE, DTEXT_CENTER, DTEXT_MIDDLE,
	          "RMNChem  --  analyse de CaHbOcNd", -1);
	dtext_opt(W / 2, 18, C_RGB(22, 22, 26), C_NONE, DTEXT_CENTER,
	          DTEXT_MIDDLE, "un indice nul supprime l'element", -1);

	for(i = 0; i < 4; i++) {
		int y = 40 + 22 * i, sel = (i == cur_field), n;
		if(sel) drect(20, y - 3, 190, y + 13, C_RGB(26, 30, 31));
		dtext(28, y, C_BLACK, lname[i]);
		n = u2s(buf, *val[i]);
		buf[n] = 0;
		dtext_opt(150, y + 4, sel ? C_BLUE : C_BLACK, C_NONE,
		          DTEXT_RIGHT, DTEXT_MIDDLE, buf, -1);
		if(sel) dtext(168, y, C_BLUE, "<");
	}

	fmt_formula(buf, in_a, in_b, in_c, in_d);
	dtext(212, 36, C_GREY, "Formule brute");
	dtext(212, 50, C_BLACK, buf);

	st = compute_di(in_a, in_b, in_c, in_d);
	dtext(212, 72, C_GREY, "Degre d'insaturation");
	if(st == DI_OK) {
		int n = u2s(buf, di_x2 / 2);
		buf[n] = 0;
		dtext(212, 86, C_BLUE, buf);
		search_db(in_a, in_b, in_c, in_d);
		n = u2s(buf, nhits);
		buf[n] = 0;
		dtext(212, 108, C_GREY, "Isomeres tabules");
		dtext(212, 122, C_BLACK, buf);
	} else if(st == DI_TOOMANY_H) {
		dtext(212, 86, C_RED, "trop de H");
	} else if(st == DI_HALF) {
		dtext(212, 86, C_RED, "demi-entier");
	} else if(st == DI_NEG) {
		dtext(212, 86, C_RED, "negatif");
	} else {
		dtext(212, 86, C_GREY, "--");
	}

	dtext(20, 146, C_GREY, "HAUT / BAS : changer de champ");
	dtext(20, 158, C_GREY, "GAUCHE / DROITE : -1 / +1    0-9 : saisie");
	dtext(20, 170, C_GREY, "DEL : remise a zero     EXE : analyser");

	draw_fkeys("C2H6O", "C4H8O2", "C7H8", "C2H5NO", "Raz", "Zones");
	dupdate();
}

/* ================================================================== */
/* Moteur de secours : formules absentes de la base                    */
/* ================================================================== */
static void screen_fallback(int a, int b, int c, int d, int st)
{
	char buf[32];
	int di = di_x2 / 2, y = MID_Y0, n;

	dclear(C_WHITE);
	draw_header(a, b, c, d,
	            st == DI_OK ? "non tabulee" : "formule invalide",
	            st == DI_OK ? "estimation par fragments"
	                        : "verifier les indices saisis");

	if(st != DI_OK) {
		dtext(10, 50, C_RED, "Formule impossible pour une molecule");
		dtext(10, 62, C_RED, "neutre a couches completes.");
		if(st == DI_TOOMANY_H) {
			dtext(10, 84, C_BLACK, "b depasse 2a + 2 + d : le squelette");
			dtext(10, 96, C_BLACK, "ne peut pas porter autant de H.");
		} else if(st == DI_HALF) {
			dtext(10, 84, C_BLACK, "DI demi-entier : b et d sont de");
			dtext(10, 96, C_BLACK, "parites incompatibles. Il s'agit");
			dtext(10, 108, C_BLACK, "d'un ion ou d'un radical.");
		} else {
			dtext(10, 84, C_BLACK, "Saisir au moins un atome.");
		}
		draw_fkeys("Retour", 0, 0, 0, 0, 0);
		dupdate();
		return;
	}

	dtext(8, y, C_BLUE, "Degre d'insaturation");
	n = u2s(buf, di);
	buf[n] = 0;
	dtext(8, y + 14, C_BLACK, "DI =");
	dtext(48, y + 14, C_BLUE, buf);
	dtext(8, y + 26, C_GREY, "DI = a+1+(d-b)/2");
	y += 44;

	dtext(8, y, C_BLUE, "Fragments compatibles");
	y += 13;

	if(di == 0)
		{ dtext(10, y, C_BLACK, "Chaine saturee, sans"); y += 11;
		  dtext(10, y, C_BLACK, "cycle ni insaturation"); y += 11; }
	if(di >= 4 && a >= 6)
		{ dtext(10, y, C_BLACK, "Noyau benzenique (DI 4)"); y += 11; }
	if(di >= 1 && di <= 3)
		{ dtext(10, y, C_BLACK, "C=C, C=O ou cycle"); y += 11; }
	if(c == 1 && di >= 1)
		{ dtext(10, y, C_BLACK, "C=O aldehyde ou cetone"); y += 11; }
	if(c == 1 && di == 0)
		{ dtext(10, y, C_BLACK, "Alcool ou ether-oxyde"); y += 11; }
	if(c >= 2 && di >= 1)
		{ dtext(10, y, C_BLACK, "Acide ou ester"); y += 11; }
	if(c >= 2 && di == 0)
		{ dtext(10, y, C_BLACK, "Diol ou diether"); y += 11; }
	if(d >= 1 && c == 0 && di == 0)
		{ dtext(10, y, C_BLACK, "Amine"); y += 11; }
	if(d >= 1 && c >= 1 && di >= 1)
		{ dtext(10, y, C_BLACK, "Amide ou derive nitre"); y += 11; }
	if(d >= 1 && c == 0 && di >= 2)
		{ dtext(10, y, C_BLACK, "Nitrile : C#N vaut 2"); y += 11; }
	if(d >= 1 && c >= 2 && di >= 1 && b >= 5)
		{ dtext(10, y, C_BLACK, "Acide amine plausible"); y += 11; }

	dtext(FB_X, MID_Y0, C_GREY, "Pas de structure");
	dtext(FB_X, MID_Y0 + 12, C_GREY, "tabulee pour cette");
	dtext(FB_X, MID_Y0 + 24, C_GREY, "formule. L'echelle");
	dtext(FB_X, MID_Y0 + 36, C_GREY, "ci-dessous rappelle");
	dtext(FB_X, MID_Y0 + 48, C_GREY, "les zones utiles.");
	n = u2s(buf, b);
	buf[n] = 0;
	dtext(FB_X, MID_Y0 + 68, C_BLACK, "H a placer :");
	dtext(FB_X + 104, MID_Y0 + 68, C_BLUE, buf);

	dline(FB_X - 6, MID_Y0 - 2, FB_X - 6, MID_Y1, C_GREY);
	dline(0, MID_Y1 + 2, W - 1, MID_Y1 + 2, C_GREY);
	draw_zone_axis();
	draw_fkeys("Retour", 0, 0, 0, 0, "Zones");
	dupdate();
}

/* ================================================================== */
/* Ecran de resultat                                                   */
/* ================================================================== */
static void screen_result(void)
{
	const molecule_t *m = &DB[hits[(int)cur_hit]];
	char buf[32];
	int n;

	dclear(C_WHITE);
	draw_header(m->a, m->b, m->c, m->d, m->name, m->family);

	n = u2s(buf, di_x2 / 2);
	buf[n] = 0;
	dtext(W - 104, 2, C_WHITE, "DI =");
	dtext(W - 64, 2, C_RGB(31, 24, 0), buf);
	if(nhits > 1) {
		int k = 0;
		buf[k++] = 'i'; buf[k++] = 's'; buf[k++] = 'o';
		buf[k++] = '.'; buf[k++] = ' ';
		k += u2s(buf + k, cur_hit + 1);
		buf[k++] = '/';
		k += u2s(buf + k, nhits);
		buf[k] = 0;
		dtext_opt(W - 6, 17, C_RGB(28, 28, 20), C_NONE,
		          DTEXT_RIGHT, DTEXT_MIDDLE, buf, -1);
	}

	dline(TABLE_X0 - 2, MID_Y0 - 2, TABLE_X0 - 2, MID_Y1, C_GREY);
	dline(0, MID_Y1 + 2, W - 1, MID_Y1 + 2, C_GREY);

	draw_molecule(m, STRUCT_X0, MID_Y0, STRUCT_X1, MID_Y1 - 4);
	draw_table(m->peaks, m->npeaks);
	draw_spectrum(m->peaks, m->npeaks);

	draw_fkeys(nhits > 1 ? "<Isom" : 0, nhits > 1 ? "Isom>" : 0,
	           show_integral ? "Integ*" : "Integ", 0, "Zones", "Retour");
	dupdate();
}

/* ================================================================== */
/* Aide-memoire des zones                                              */
/* ================================================================== */
static void screen_zones(void)
{
	dclear(C_WHITE);
	drect(0, 0, W - 1, 20, C_RGB(0, 0, 14));
	dtext(6, 5, C_WHITE, "Deplacements chimiques -- RMN du proton");

	dtext(8, 28, C_BLUE, "Champ fort");
	dtext(8,  42, C_BLACK, "0,8-1,3  CH3 alcane");
	dtext(8,  54, C_BLACK, "1,2-1,6  CH2 alcane");
	dtext(8,  66, C_BLACK, "1,4-1,8  CH alcane");
	dtext(8,  78, C_BLACK, "2,0-2,6  CH en a d'un C=O");
	dtext(8,  90, C_BLACK, "2,2-3,0  CH-Ar");
	dtext(8, 102, C_BLACK, "2,2-2,9  CH-N amine");
	dtext(8, 114, C_BLACK, "3,3-4,3  CH-O");

	dtext(204, 28, C_BLUE, "Champ faible");
	dtext(204,  42, C_BLACK, "4,5-6,5  =CH vinylique");
	dtext(204,  54, C_BLACK, "6,5-8,5  H aromatique");
	dtext(204,  66, C_BLACK, "9,5-10,5 CHO aldehyde");
	dtext(204,  78, C_BLACK, "10-13    COOH acide");
	dtext(204,  90, C_BLACK, "0,5-5,5  OH alcool");
	dtext(204, 102, C_BLACK, "4-8      OH phenol");
	dtext(204, 114, C_BLACK, "1-5 NH amine, 5-8 amide");

	dline(6, 130, W - 6, 130, C_GREY);
	dtext(8, 136, C_RED, "Protons labiles : OH, NH, COOH");
	dtext(8, 148, C_BLACK, "Signal large, non couple, position tres");
	dtext(8, 160, C_BLACK, "variable. L'echange avec D2O l'efface :");
	dtext(8, 172, C_BLACK, "c'est le test qui permet de le reperer.");
	dtext(8, 186, C_RED, "Multiplicite : regle des (n+1) uplets");
	dtext(8, 198, C_BLACK, "n = H portes par les carbones voisins.");

	draw_fkeys("Retour", 0, 0, 0, 0, 0);
	dupdate();
}

/* ================================================================== */
/* Boucle principale                                                   */
/* ================================================================== */
enum { SC_INPUT, SC_RESULT, SC_FALLBACK, SC_ZONES };

static int screen = SC_INPUT, prev_screen = SC_INPUT;

static void redraw(void)
{
	switch(screen) {
	case SC_INPUT:  screen_input();  break;
	case SC_RESULT: screen_result(); break;
	case SC_ZONES:  screen_zones();  break;
	case SC_FALLBACK:
		screen_fallback(in_a, in_b, in_c, in_d,
		                compute_di(in_a, in_b, in_c, in_d));
		break;
	}
}

int main(void)
{
	uint8_t *val[4] = { &in_a, &in_b, &in_c, &in_d };
	static const uint8_t vmax[4] = { 40, 80, 12, 12 };

	while(1) {
		key_event_t ev;
		int k;

		redraw();
		ev = getkey();
		k = ev.key;

		if(k == KEY_EXIT) {
			if(screen == SC_INPUT) break;
			screen = (screen == SC_ZONES) ? prev_screen : SC_INPUT;
			continue;
		}
		if(screen == SC_ZONES) {
			if(k == KEY_F1) screen = prev_screen;
			continue;
		}
		if(screen == SC_FALLBACK) {
			if(k == KEY_F1) screen = SC_INPUT;
			if(k == KEY_F6) { prev_screen = SC_FALLBACK; screen = SC_ZONES; }
			continue;
		}
		if(screen == SC_RESULT && k == KEY_F5) {
			prev_screen = SC_RESULT; screen = SC_ZONES; continue;
		}
		if(screen == SC_INPUT && k == KEY_F6) {
			prev_screen = SC_INPUT; screen = SC_ZONES; continue;
		}

		if(screen == SC_INPUT) {
			if(k == KEY_UP)
				{ cur_field = (int8_t)((cur_field + 3) % 4); typing = 0; }
			if(k == KEY_DOWN)
				{ cur_field = (int8_t)((cur_field + 1) % 4); typing = 0; }
			if(k == KEY_LEFT && *val[(int)cur_field] > 0)
				{ (*val[(int)cur_field])--; typing = 0; }
			if(k == KEY_RIGHT && *val[(int)cur_field] < vmax[(int)cur_field])
				{ (*val[(int)cur_field])++; typing = 0; }
			if(k == KEY_DEL) { *val[(int)cur_field] = 0; typing = 0; }

			if(k >= KEY_0 && k <= KEY_9) {
				int digit = k - KEY_0;
				int v = typing ? *val[(int)cur_field] * 10 + digit : digit;
				if(v > vmax[(int)cur_field]) v = digit;
				*val[(int)cur_field] = (uint8_t)v;
				typing = 1;
			}

			if(k == KEY_F1) { in_a=2; in_b=6; in_c=1; in_d=0; typing=0; }
			if(k == KEY_F2) { in_a=4; in_b=8; in_c=2; in_d=0; typing=0; }
			if(k == KEY_F3) { in_a=7; in_b=8; in_c=0; in_d=0; typing=0; }
			if(k == KEY_F4) { in_a=2; in_b=5; in_c=1; in_d=1; typing=0; }
			if(k == KEY_F5) { in_a=in_b=in_c=in_d=0; typing=0; }

			if(k == KEY_EXE) {
				int st = compute_di(in_a, in_b, in_c, in_d);
				typing = 0;
				search_db(in_a, in_b, in_c, in_d);
				screen = (st == DI_OK && nhits > 0) ? SC_RESULT
				                                    : SC_FALLBACK;
			}
			continue;
		}

		if(screen == SC_RESULT) {
			if((k == KEY_LEFT || k == KEY_F1) && nhits > 1)
				cur_hit = (int8_t)((cur_hit + nhits - 1) % nhits);
			if((k == KEY_RIGHT || k == KEY_F2) && nhits > 1)
				cur_hit = (int8_t)((cur_hit + 1) % nhits);
			if(k == KEY_F3) show_integral = !show_integral;
			if(k == KEY_F6) screen = SC_INPUT;
		}
	}

	return 1;
}
