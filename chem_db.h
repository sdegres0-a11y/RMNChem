/*
 * chem_db.h -- Base de donnees chimique compacte pour RMNChem
 *
 * Toutes les donnees sont "const" : elles restent en ROM (section .rodata du
 * .g3a) et ne consomment pas de RAM sur la calculatrice.
 *
 * Conventions de compacite :
 *   - deplacements chimiques stockes en ppm x 10 sur 1 octet (0.0 -> 25.5 ppm)
 *   - coordonnees 2D en pixels relatifs sur 1 octet signe (-128 .. +127)
 *   - element / nombre d'hydrogenes / drapeaux d'affichage en champs de bits
 *   - attributions stockees comme index dans une table de chaines partagee
 *
 * Cout memoire : 3 octets / atome, 2 octets / liaison, 3 octets / signal RMN.
 */

#ifndef CHEM_DB_H
#define CHEM_DB_H

#include <stdint.h>

/* ------------------------------------------------------------------ */
/* Elements                                                            */
/* ------------------------------------------------------------------ */
enum {
	EL_C = 0,
	EL_O = 1,
	EL_N = 2
};

/* ------------------------------------------------------------------ */
/* Multiplicites                                                       */
/* ------------------------------------------------------------------ */
enum {
	M_S = 0,	/* singulet     */
	M_D = 1,	/* doublet      */
	M_T = 2,	/* triplet      */
	M_Q = 3,	/* quadruplet   */
	M_QT = 4,	/* quintuplet   */
	M_SX = 5,	/* sextuplet    */
	M_SP = 6,	/* septuplet    */
	M_M = 7		/* multiplet    */
};

/* Nombre de raies dessinees pour chaque multiplicite */
static const uint8_t MULT_LINES[8] = { 1, 2, 3, 4, 5, 6, 7, 3 };
static const char * const MULT_TXT[8] = { "s", "d", "t", "q", "quint",
                                          "sext", "sept", "m" };

/* ------------------------------------------------------------------ */
/* Atome : 3 octets                                                    */
/* ------------------------------------------------------------------ */
typedef struct {
	int8_t x, y;		/* position relative, en pixels                */
	uint8_t el   : 2;	/* EL_C / EL_O / EL_N                          */
	uint8_t nH   : 3;	/* hydrogenes portes (0..4)                    */
	uint8_t show : 1;	/* 1 = ecrire le label, 0 = sommet implicite   */
	uint8_t hl   : 1;	/* 1 = groupe caracteristique (trace en rouge) */
	uint8_t pad  : 1;
} atom_t;

/* ------------------------------------------------------------------ */
/* Liaison : 2 octets                                                  */
/* ------------------------------------------------------------------ */
typedef struct {
	uint8_t i : 4;		/* index du premier atome  (max 16 atomes)     */
	uint8_t j : 4;		/* index du second atome                       */
	uint8_t order : 2;	/* 1 simple, 2 double, 3 triple                */
	uint8_t arom  : 1;	/* liaison d'un cycle aromatique               */
	uint8_t pad   : 5;
} bond_t;

/* ------------------------------------------------------------------ */
/* Signal RMN : 3 octets                                               */
/* ------------------------------------------------------------------ */
typedef struct {
	uint8_t shift;		/* delta x 10 : 72 => 7,2 ppm                  */
	uint8_t nH   : 4;	/* integration (1..15 H)                       */
	uint8_t mult : 3;	/* M_S .. M_M                                  */
	uint8_t exch : 1;	/* proton labile : signal large, disparait D2O */
	uint8_t attr;		/* index dans ATTR[]                           */
} nmr_peak_t;

/* ------------------------------------------------------------------ */
/* Molecule                                                            */
/* ------------------------------------------------------------------ */
typedef struct {
	uint8_t a, b, c, d;		/* indices de CaHbOcNd                 */
	const char *name;		/* nom d'usage (celui attendu en DS)   */
	const char *family;		/* famille / remarque courte           */
	const atom_t *atoms;
	const bond_t *bonds;
	const nmr_peak_t *peaks;
	uint8_t natoms, nbonds, npeaks;
} molecule_t;

/* ------------------------------------------------------------------ */
/* Table d'attributions partagee                                       */
/* ------------------------------------------------------------------ */
enum {
	A_NONE = 0, A_CH3, A_CH2, A_CH, A_OH, A_NH, A_NH2, A_CHO, A_COOH,
	A_AROM, A_OCH3, A_OCH2, A_CH3CO, A_CH2CO, A_CHO_LNK, A_ORTHO,
	A_META, A_PARA, A_CH2EQ, A_CHEQ, A_ALCYNE, A_NCH3, A_NCH2, A_NCH,
	A_HCOO, A_TBU, A_LABILES, A_CYCLE, A_ARCH3, A_ARCH2, A_CH4, A_CH3O
};

static const char * const ATTR[] = {
	"",              "CH3",           "CH2",          "CH",
	"OH",            "NH",            "NH2",          "CHO",
	"COOH",          "H aromatiques", "OCH3",         "OCH2",
	"CH3 en a C=O",  "CH2 en a C=O",  "CH-O",         "H ortho",
	"H meta",        "H para",        "=CH2",         "=CH-",
	"H alcyne",      "CH3-N",         "CH2-N",        "CH-N",
	"H du formiate", "C(CH3)3",       "H labiles",    "H du cycle",
	"CH3-Ar",        "CH2-Ar",        "CH4",          "CH3-O"
};

/* ================================================================== */
/* Aides de mise en page                                               */
/*                                                                     */
/* Chaine en zigzag : longueur de liaison ~22 px (dx = 19, dy = 11).   */
/* Hexagone aromatique : rayon 22 px, sommet 0 a droite.               */
/* Le rendu recentre automatiquement, seules les positions relatives   */
/* comptent.                                                           */
/* ================================================================== */
#define ZX(i) ((int8_t)(19 * (i)))
#define ZY(i) ((int8_t)(((i) & 1) ? 0 : 11))

/* Sommets de l'hexagone, dans l'ordre 0..5 */
#define HEX6 \
	{  22,   0, EL_C, 1, 0, 0, 0 }, \
	{  11, -19, EL_C, 1, 0, 0, 0 }, \
	{ -11, -19, EL_C, 1, 0, 0, 0 }, \
	{ -22,   0, EL_C, 1, 0, 0, 0 }, \
	{ -11,  19, EL_C, 1, 0, 0, 0 }, \
	{  11,  19, EL_C, 1, 0, 0, 0 }

/* Idem, mais le sommet 0 porte le substituant (donc 0 H) */
#define HEX5 \
	{  22,   0, EL_C, 0, 0, 0, 0 }, \
	{  11, -19, EL_C, 1, 0, 0, 0 }, \
	{ -11, -19, EL_C, 1, 0, 0, 0 }, \
	{ -22,   0, EL_C, 1, 0, 0, 0 }, \
	{ -11,  19, EL_C, 1, 0, 0, 0 }, \
	{  11,  19, EL_C, 1, 0, 0, 0 }

/* Cycle disubstitue en para : sommets 0 et 3 nus */
#define HEX4P \
	{  22,   0, EL_C, 0, 0, 0, 0 }, \
	{  11, -19, EL_C, 1, 0, 0, 0 }, \
	{ -11, -19, EL_C, 1, 0, 0, 0 }, \
	{ -22,   0, EL_C, 0, 0, 0, 0 }, \
	{ -11,  19, EL_C, 1, 0, 0, 0 }, \
	{  11,  19, EL_C, 1, 0, 0, 0 }

/* Cycle disubstitue en ortho : sommets 0 et 1 nus */
#define HEX4O \
	{  22,   0, EL_C, 0, 0, 0, 0 }, \
	{  11, -19, EL_C, 0, 0, 0, 0 }, \
	{ -11, -19, EL_C, 1, 0, 0, 0 }, \
	{ -22,   0, EL_C, 1, 0, 0, 0 }, \
	{ -11,  19, EL_C, 1, 0, 0, 0 }, \
	{  11,  19, EL_C, 1, 0, 0, 0 }

/* Les six liaisons du cycle */
#define HEXB \
	{ 0, 1, 1, 1, 0 }, { 1, 2, 1, 1, 0 }, { 2, 3, 1, 1, 0 }, \
	{ 3, 4, 1, 1, 0 }, { 4, 5, 1, 1, 0 }, { 5, 0, 1, 1, 0 }

/* ================================================================== */
/* 1. Hydrocarbures                                                    */
/* ================================================================== */

static const atom_t at_methane[] = { { 0, 0, EL_C, 4, 1, 0, 0 } };
static const nmr_peak_t pk_methane[] = { { 2, 4, M_S, 0, A_CH4 } };

static const atom_t at_ethane[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 }, { ZX(1), ZY(1), EL_C, 3, 1, 0, 0 }
};
static const bond_t bd_chain1[] = { { 0, 1, 1, 0, 0 } };
static const nmr_peak_t pk_ethane[] = { { 9, 6, M_S, 0, A_CH3 } };

static const atom_t at_ethene[] = {
	{ 0, 0, EL_C, 2, 1, 0, 1 }, { 22, 0, EL_C, 2, 1, 0, 1 }
};
static const bond_t bd_double1[] = { { 0, 1, 2, 0, 0 } };
static const nmr_peak_t pk_ethene[] = { { 53, 4, M_S, 0, A_CH2EQ } };

static const atom_t at_ethyne[] = {
	{ 0, 0, EL_C, 1, 1, 0, 1 }, { 24, 0, EL_C, 1, 1, 0, 1 }
};
static const bond_t bd_triple1[] = { { 0, 1, 3, 0, 0 } };
static const nmr_peak_t pk_ethyne[] = { { 20, 2, M_S, 0, A_ALCYNE } };

static const atom_t at_cyclohexane[] = {
	{  22,   0, EL_C, 2, 0, 0, 0 }, {  11, -19, EL_C, 2, 0, 0, 0 },
	{ -11, -19, EL_C, 2, 0, 0, 0 }, { -22,   0, EL_C, 2, 0, 0, 0 },
	{ -11,  19, EL_C, 2, 0, 0, 0 }, {  11,  19, EL_C, 2, 0, 0, 0 }
};
static const bond_t bd_ring6[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 1, 0, 0 }, { 2, 3, 1, 0, 0 },
	{ 3, 4, 1, 0, 0 }, { 4, 5, 1, 0, 0 }, { 5, 0, 1, 0, 0 }
};
static const nmr_peak_t pk_cyclohexane[] = { { 14, 12, M_S, 0, A_CYCLE } };

static const atom_t at_benzene[] = { HEX6 };
static const bond_t bd_benzene[] = { HEXB };
static const nmr_peak_t pk_benzene[] = { { 73, 6, M_S, 0, A_AROM } };

static const atom_t at_toluene[] = { HEX5, { 44, 0, EL_C, 3, 1, 0, 0 } };
static const bond_t bd_sub1[] = { HEXB, { 0, 6, 1, 0, 0 } };
static const nmr_peak_t pk_toluene[] = {
	{ 24, 3, M_S, 0, A_ARCH3 }, { 72, 5, M_M, 0, A_AROM }
};

static const atom_t at_ethylbenzene[] = {
	HEX5, { 44, 0, EL_C, 2, 1, 0, 0 }, { 63, -11, EL_C, 3, 1, 0, 0 }
};
static const bond_t bd_sub2[] = { HEXB, { 0, 6, 1, 0, 0 }, { 6, 7, 1, 0, 0 } };
static const nmr_peak_t pk_ethylbenzene[] = {
	{ 13, 3, M_T, 0, A_CH3 }, { 27, 2, M_Q, 0, A_ARCH2 },
	{ 72, 5, M_M, 0, A_AROM }
};

static const atom_t at_pxylene[] = {
	HEX4P, { 44, 0, EL_C, 3, 1, 0, 0 }, { -44, 0, EL_C, 3, 1, 0, 0 }
};
static const bond_t bd_pxylene[] = {
	HEXB, { 0, 6, 1, 0, 0 }, { 3, 7, 1, 0, 0 }
};
static const nmr_peak_t pk_pxylene[] = {
	{ 23, 6, M_S, 0, A_ARCH3 }, { 71, 4, M_S, 0, A_AROM }
};

static const atom_t at_oxylene[] = {
	HEX4O, { 44, 0, EL_C, 3, 1, 0, 0 }, { 22, -38, EL_C, 3, 1, 0, 0 }
};
static const bond_t bd_oxylene[] = {
	HEXB, { 0, 6, 1, 0, 0 }, { 1, 7, 1, 0, 0 }
};
static const nmr_peak_t pk_oxylene[] = {
	{ 23, 6, M_S, 0, A_ARCH3 }, { 71, 4, M_M, 0, A_AROM }
};

static const atom_t at_styrene[] = {
	HEX5, { 44, 0, EL_C, 1, 1, 0, 1 }, { 63, -11, EL_C, 2, 1, 0, 1 }
};
static const bond_t bd_styrene[] = {
	HEXB, { 0, 6, 1, 0, 0 }, { 6, 7, 2, 0, 0 }
};
static const nmr_peak_t pk_styrene[] = {
	{ 53, 1, M_D, 0, A_CH2EQ }, { 58, 1, M_D, 0, A_CH2EQ },
	{ 67, 1, M_M, 0, A_CHEQ  }, { 73, 5, M_M, 0, A_AROM }
};

/* ================================================================== */
/* 2. Alcools et ethers                                                */
/* ================================================================== */

static const atom_t at_methanol[] = {
	{ 0, 0, EL_C, 3, 1, 0, 0 }, { 24, 0, EL_O, 1, 1, 1, 0 }
};
static const nmr_peak_t pk_methanol[] = {
	{ 34, 3, M_S, 0, A_CH3O }, { 20, 1, M_S, 1, A_OH }
};

static const atom_t at_ethanol[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_C, 2, 1, 0, 0 },
	{ ZX(2), ZY(2), EL_O, 1, 1, 1, 0 }
};
static const bond_t bd_chain2[] = { { 0, 1, 1, 0, 0 }, { 1, 2, 1, 0, 0 } };
static const nmr_peak_t pk_ethanol[] = {
	{ 12, 3, M_T, 0, A_CH3 }, { 23, 1, M_S, 1, A_OH },
	{ 37, 2, M_Q, 0, A_OCH2 }
};

static const atom_t at_dimethylether[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_O, 0, 1, 1, 0 },
	{ ZX(2), ZY(2), EL_C, 3, 1, 0, 0 }
};
static const nmr_peak_t pk_dimethylether[] = { { 33, 6, M_S, 0, A_OCH3 } };

static const atom_t at_propan1ol[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_C, 2, 0, 0, 0 },
	{ ZX(2), ZY(2), EL_C, 2, 0, 0, 0 },
	{ ZX(3), ZY(3), EL_O, 1, 1, 1, 0 }
};
static const bond_t bd_chain3[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 1, 0, 0 }, { 2, 3, 1, 0, 0 }
};
static const nmr_peak_t pk_propan1ol[] = {
	{  9, 3, M_T,  0, A_CH3 }, { 16, 2, M_SX, 0, A_CH2 },
	{ 25, 1, M_S,  1, A_OH  }, { 36, 2, M_T,  0, A_OCH2 }
};

static const atom_t at_propan2ol[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_C, 1, 0, 0, 0 },
	{ ZX(2), ZY(2), EL_C, 3, 1, 0, 0 },
	{ 19, -22, EL_O, 1, 1, 1, 0 }
};
static const bond_t bd_iso3[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 1, 0, 0 }, { 1, 3, 1, 0, 0 }
};
static const nmr_peak_t pk_propan2ol[] = {
	{ 12, 6, M_D, 0, A_CH3 }, { 22, 1, M_S, 1, A_OH },
	{ 39, 1, M_SP, 0, A_CH  }
};

static const atom_t at_butan1ol[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_C, 2, 0, 0, 0 },
	{ ZX(2), ZY(2), EL_C, 2, 0, 0, 0 },
	{ ZX(3), ZY(3), EL_C, 2, 0, 0, 0 },
	{ ZX(4), ZY(4), EL_O, 1, 1, 1, 0 }
};
static const bond_t bd_chain4[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 1, 0, 0 },
	{ 2, 3, 1, 0, 0 }, { 3, 4, 1, 0, 0 }
};
static const nmr_peak_t pk_butan1ol[] = {
	{  9, 3, M_T,  0, A_CH3 }, { 14, 2, M_SX, 0, A_CH2 },
	{ 16, 2, M_QT, 0, A_CH2 }, { 20, 1, M_S,  1, A_OH  },
	{ 36, 2, M_T,  0, A_OCH2 }
};

static const atom_t at_tbutanol[] = {
	{   0,   0, EL_C, 0, 0, 0, 0 },
	{ -22,   0, EL_C, 3, 1, 0, 0 },
	{   0, -22, EL_C, 3, 1, 0, 0 },
	{   0,  22, EL_C, 3, 1, 0, 0 },
	{  24,   0, EL_O, 1, 1, 1, 0 }
};
static const bond_t bd_quat4[] = {
	{ 0, 1, 1, 0, 0 }, { 0, 2, 1, 0, 0 },
	{ 0, 3, 1, 0, 0 }, { 0, 4, 1, 0, 0 }
};
static const nmr_peak_t pk_tbutanol[] = {
	{ 13, 9, M_S, 0, A_TBU }, { 20, 1, M_S, 1, A_OH }
};

static const atom_t at_diethylether[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_C, 2, 0, 0, 0 },
	{ ZX(2), ZY(2), EL_O, 0, 1, 1, 0 },
	{ ZX(3), ZY(3), EL_C, 2, 0, 0, 0 },
	{ ZX(4), ZY(4), EL_C, 3, 1, 0, 0 }
};
static const nmr_peak_t pk_diethylether[] = {
	{ 12, 6, M_T, 0, A_CH3 }, { 34, 4, M_Q, 0, A_OCH2 }
};

static const atom_t at_phenol[] = { HEX5, { 46, 0, EL_O, 1, 1, 1, 0 } };
static const nmr_peak_t pk_phenol[] = {
	{ 54, 1, M_S, 1, A_OH   }, { 68, 3, M_M, 0, A_ORTHO },
	{ 72, 2, M_M, 0, A_META }
};

static const atom_t at_anisole[] = {
	HEX5, { 46, 0, EL_O, 0, 1, 1, 0 }, { 65, -11, EL_C, 3, 1, 0, 0 }
};
static const nmr_peak_t pk_anisole[] = {
	{ 38, 3, M_S, 0, A_OCH3 }, { 69, 3, M_M, 0, A_ORTHO },
	{ 73, 2, M_M, 0, A_META }
};

/* ================================================================== */
/* 3. Aldehydes et cetones                                             */
/* ================================================================== */

static const atom_t at_ethanal[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_C, 1, 0, 0, 1 },
	{ ZX(2), ZY(2), EL_O, 0, 1, 1, 0 }
};
static const bond_t bd_carbonyl_end2[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 2, 0, 0 }
};
static const nmr_peak_t pk_ethanal[] = {
	{ 22, 3, M_D, 0, A_CH3CO }, { 98, 1, M_Q, 0, A_CHO }
};

static const atom_t at_propanal[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_C, 2, 0, 0, 0 },
	{ ZX(2), ZY(2), EL_C, 1, 0, 0, 1 },
	{ ZX(3), ZY(3), EL_O, 0, 1, 1, 0 }
};
static const bond_t bd_carbonyl_end3[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 1, 0, 0 }, { 2, 3, 2, 0, 0 }
};
static const nmr_peak_t pk_propanal[] = {
	{ 11, 3, M_T, 0, A_CH3 }, { 24, 2, M_Q, 0, A_CH2CO },
	{ 98, 1, M_T, 0, A_CHO }
};

static const atom_t at_butanal[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_C, 2, 0, 0, 0 },
	{ ZX(2), ZY(2), EL_C, 2, 0, 0, 0 },
	{ ZX(3), ZY(3), EL_C, 1, 0, 0, 1 },
	{ ZX(4), ZY(4), EL_O, 0, 1, 1, 0 }
};
static const bond_t bd_carbonyl_end4[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 1, 0, 0 },
	{ 2, 3, 1, 0, 0 }, { 3, 4, 2, 0, 0 }
};
static const nmr_peak_t pk_butanal[] = {
	{ 10, 3, M_T,  0, A_CH3 }, { 17, 2, M_SX, 0, A_CH2 },
	{ 24, 2, M_T,  0, A_CH2CO }, { 98, 1, M_T, 0, A_CHO }
};

static const atom_t at_propanone[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_C, 0, 0, 0, 1 },
	{ ZX(2), ZY(2), EL_C, 3, 1, 0, 0 },
	{ 19, -22, EL_O, 0, 1, 1, 0 }
};
static const bond_t bd_ketone3[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 1, 0, 0 }, { 1, 3, 2, 0, 0 }
};
static const nmr_peak_t pk_propanone[] = { { 21, 6, M_S, 0, A_CH3CO } };

static const atom_t at_butanone[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_C, 0, 0, 0, 1 },
	{ ZX(2), ZY(2), EL_C, 2, 0, 0, 0 },
	{ ZX(3), ZY(3), EL_C, 3, 1, 0, 0 },
	{ 19, -22, EL_O, 0, 1, 1, 0 }
};
static const bond_t bd_ketone4[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 1, 0, 0 },
	{ 2, 3, 1, 0, 0 }, { 1, 4, 2, 0, 0 }
};
static const nmr_peak_t pk_butanone[] = {
	{ 10, 3, M_T, 0, A_CH3 }, { 21, 3, M_S, 0, A_CH3CO },
	{ 24, 2, M_Q, 0, A_CH2CO }
};

static const atom_t at_benzaldehyde[] = {
	HEX5, { 44, 0, EL_C, 1, 0, 0, 1 }, { 63, -11, EL_O, 0, 1, 1, 0 }
};
static const bond_t bd_ar_carbonyl[] = {
	HEXB, { 0, 6, 1, 0, 0 }, { 6, 7, 2, 0, 0 }
};
static const nmr_peak_t pk_benzaldehyde[] = {
	{ 75, 3, M_M, 0, A_META }, { 79, 2, M_M, 0, A_ORTHO },
	{ 100, 1, M_S, 0, A_CHO }
};

static const atom_t at_acetophenone[] = {
	HEX5, { 44, 0, EL_C, 0, 0, 0, 1 }, { 44, -22, EL_O, 0, 1, 1, 0 },
	{ 63, 11, EL_C, 3, 1, 0, 0 }
};
static const bond_t bd_ar_ketone[] = {
	HEXB, { 0, 6, 1, 0, 0 }, { 6, 7, 2, 0, 0 }, { 6, 8, 1, 0, 0 }
};
static const nmr_peak_t pk_acetophenone[] = {
	{ 26, 3, M_S, 0, A_CH3CO }, { 75, 3, M_M, 0, A_META },
	{ 79, 2, M_M, 0, A_ORTHO }
};

/* ================================================================== */
/* 4. Acides carboxyliques et esters                                   */
/* ================================================================== */

static const atom_t at_methanoic[] = {
	{  0,   0, EL_C, 1, 0, 0, 1 },
	{ 19, -11, EL_O, 0, 1, 1, 0 },
	{ 19,  11, EL_O, 1, 1, 1, 0 }
};
static const bond_t bd_acid_head[] = { { 0, 1, 2, 0, 0 }, { 0, 2, 1, 0, 0 } };
static const nmr_peak_t pk_methanoic[] = {
	{ 81, 1, M_S, 0, A_CHO }, { 114, 1, M_S, 1, A_COOH }
};

static const atom_t at_ethanoic[] = {
	{   0,  11, EL_C, 3, 1, 0, 0 },
	{  19,   0, EL_C, 0, 0, 0, 1 },
	{  19, -22, EL_O, 0, 1, 1, 0 },
	{  38,  11, EL_O, 1, 1, 1, 0 }
};
static const bond_t bd_acid2[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 2, 0, 0 }, { 1, 3, 1, 0, 0 }
};
static const nmr_peak_t pk_ethanoic[] = {
	{ 21, 3, M_S, 0, A_CH3CO }, { 114, 1, M_S, 1, A_COOH }
};

static const atom_t at_propanoic[] = {
	{   0,  11, EL_C, 3, 1, 0, 0 },
	{  19,   0, EL_C, 2, 0, 0, 0 },
	{  38,  11, EL_C, 0, 0, 0, 1 },
	{  38, -11, EL_O, 0, 1, 1, 0 },
	{  57,  22, EL_O, 1, 1, 1, 0 }
};
static const bond_t bd_acid3[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 1, 0, 0 },
	{ 2, 3, 2, 0, 0 }, { 2, 4, 1, 0, 0 }
};
static const nmr_peak_t pk_propanoic[] = {
	{ 12, 3, M_T, 0, A_CH3 }, { 24, 2, M_Q, 0, A_CH2CO },
	{ 116, 1, M_S, 1, A_COOH }
};

static const atom_t at_benzoic[] = {
	HEX5, { 44, 0, EL_C, 0, 0, 0, 1 },
	{ 63, -11, EL_O, 0, 1, 1, 0 }, { 63, 11, EL_O, 1, 1, 1, 0 }
};
static const bond_t bd_ar_acid[] = {
	HEXB, { 0, 6, 1, 0, 0 }, { 6, 7, 2, 0, 0 }, { 6, 8, 1, 0, 0 }
};
static const nmr_peak_t pk_benzoic[] = {
	{ 75, 3, M_M, 0, A_META }, { 81, 2, M_M, 0, A_ORTHO },
	{ 120, 1, M_S, 1, A_COOH }
};

static const atom_t at_methylformate[] = {
	{   0,   0, EL_C, 1, 0, 0, 1 },
	{ -19, -11, EL_O, 0, 1, 1, 0 },
	{  19,  11, EL_O, 0, 1, 1, 0 },
	{  38,   0, EL_C, 3, 1, 0, 0 }
};
static const bond_t bd_formate[] = {
	{ 0, 1, 2, 0, 0 }, { 0, 2, 1, 0, 0 }, { 2, 3, 1, 0, 0 }
};
static const nmr_peak_t pk_methylformate[] = {
	{ 38, 3, M_S, 0, A_OCH3 }, { 81, 1, M_S, 0, A_HCOO }
};

static const atom_t at_methylacetate[] = {
	{   0,  11, EL_C, 3, 1, 0, 0 },
	{  19,   0, EL_C, 0, 0, 0, 1 },
	{  19, -22, EL_O, 0, 1, 1, 0 },
	{  38,  11, EL_O, 0, 1, 1, 0 },
	{  57,   0, EL_C, 3, 1, 0, 0 }
};
static const bond_t bd_ester4[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 2, 0, 0 },
	{ 1, 3, 1, 0, 0 }, { 3, 4, 1, 0, 0 }
};
static const nmr_peak_t pk_methylacetate[] = {
	{ 21, 3, M_S, 0, A_CH3CO }, { 37, 3, M_S, 0, A_OCH3 }
};

static const atom_t at_ethylacetate[] = {
	{   0,  11, EL_C, 3, 1, 0, 0 },
	{  19,   0, EL_C, 0, 0, 0, 1 },
	{  19, -22, EL_O, 0, 1, 1, 0 },
	{  38,  11, EL_O, 0, 1, 1, 0 },
	{  57,   0, EL_C, 2, 0, 0, 0 },
	{  76,  11, EL_C, 3, 1, 0, 0 }
};
static const bond_t bd_ester5[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 2, 0, 0 }, { 1, 3, 1, 0, 0 },
	{ 3, 4, 1, 0, 0 }, { 4, 5, 1, 0, 0 }
};
static const nmr_peak_t pk_ethylacetate[] = {
	{ 13, 3, M_T, 0, A_CH3 }, { 21, 3, M_S, 0, A_CH3CO },
	{ 41, 2, M_Q, 0, A_OCH2 }
};

static const atom_t at_methylbenzoate[] = {
	HEX5, { 44, 0, EL_C, 0, 0, 0, 1 }, { 44, -22, EL_O, 0, 1, 1, 0 },
	{ 63, 11, EL_O, 0, 1, 1, 0 }, { 82, 0, EL_C, 3, 1, 0, 0 }
};
static const bond_t bd_ar_ester[] = {
	HEXB, { 0, 6, 1, 0, 0 }, { 6, 7, 2, 0, 0 },
	{ 6, 8, 1, 0, 0 }, { 8, 9, 1, 0, 0 }
};
static const nmr_peak_t pk_methylbenzoate[] = {
	{ 39, 3, M_S, 0, A_OCH3 }, { 74, 3, M_M, 0, A_META },
	{ 80, 2, M_M, 0, A_ORTHO }
};

/* ================================================================== */
/* 5. Composes azotes                                                  */
/* ================================================================== */

static const atom_t at_methylamine[] = {
	{ 0, 0, EL_C, 3, 1, 0, 0 }, { 24, 0, EL_N, 2, 1, 1, 0 }
};
static const nmr_peak_t pk_methylamine[] = {
	{ 11, 2, M_S, 1, A_NH2 }, { 24, 3, M_S, 0, A_NCH3 }
};

static const atom_t at_ethylamine[] = {
	{ ZX(0), ZY(0), EL_C, 3, 1, 0, 0 },
	{ ZX(1), ZY(1), EL_C, 2, 0, 0, 0 },
	{ ZX(2), ZY(2), EL_N, 2, 1, 1, 0 }
};
static const nmr_peak_t pk_ethylamine[] = {
	{ 11, 3, M_T, 0, A_CH3 }, { 12, 2, M_S, 1, A_NH2 },
	{ 27, 2, M_Q, 0, A_NCH2 }
};

static const atom_t at_aniline[] = { HEX5, { 46, 0, EL_N, 2, 1, 1, 0 } };
static const nmr_peak_t pk_aniline[] = {
	{ 36, 2, M_S, 1, A_NH2 }, { 66, 2, M_M, 0, A_ORTHO },
	{ 67, 1, M_M, 0, A_PARA }, { 71, 2, M_M, 0, A_META }
};

static const atom_t at_acetonitrile[] = {
	{  0, 0, EL_C, 3, 1, 0, 0 },
	{ 24, 0, EL_C, 0, 0, 0, 1 },
	{ 48, 0, EL_N, 0, 1, 1, 0 }
};
static const bond_t bd_nitrile[] = { { 0, 1, 1, 0, 0 }, { 1, 2, 3, 0, 0 } };
static const nmr_peak_t pk_acetonitrile[] = { { 20, 3, M_S, 0, A_CH3 } };

static const atom_t at_formamide[] = {
	{  0,   0, EL_C, 1, 0, 0, 1 },
	{ 19, -11, EL_O, 0, 1, 1, 0 },
	{ 19,  11, EL_N, 2, 1, 1, 0 }
};
static const nmr_peak_t pk_formamide[] = {
	{ 65, 2, M_S, 1, A_NH2 }, { 81, 1, M_S, 0, A_CHO }
};

static const atom_t at_acetamide[] = {
	{   0,  11, EL_C, 3, 1, 0, 0 },
	{  19,   0, EL_C, 0, 0, 0, 1 },
	{  19, -22, EL_O, 0, 1, 1, 0 },
	{  38,  11, EL_N, 2, 1, 1, 0 }
};
static const nmr_peak_t pk_acetamide[] = {
	{ 20, 3, M_S, 0, A_CH3CO }, { 59, 2, M_S, 1, A_NH2 }
};

static const atom_t at_nmethylformamide[] = {
	{   0,   0, EL_C, 1, 0, 0, 1 },
	{ -19, -11, EL_O, 0, 1, 1, 0 },
	{  19,  11, EL_N, 1, 1, 1, 0 },
	{  38,   0, EL_C, 3, 1, 0, 0 }
};
static const nmr_peak_t pk_nmethylformamide[] = {
	{ 29, 3, M_D, 0, A_NCH3 }, { 64, 1, M_S, 1, A_NH },
	{ 81, 1, M_S, 0, A_CHO }
};

static const atom_t at_nitrobenzene[] = {
	HEX5, { 46, 0, EL_N, 0, 1, 1, 0 },
	{ 65, -11, EL_O, 0, 1, 1, 0 }, { 65, 11, EL_O, 0, 1, 1, 0 }
};
static const bond_t bd_nitro[] = {
	HEXB, { 0, 6, 1, 0, 0 }, { 6, 7, 2, 0, 0 }, { 6, 8, 1, 0, 0 }
};
static const nmr_peak_t pk_nitrobenzene[] = {
	{ 75, 2, M_M, 0, A_META }, { 77, 1, M_M, 0, A_PARA },
	{ 82, 2, M_M, 0, A_ORTHO }
};

static const atom_t at_glycine[] = {
	{   0,  11, EL_N, 2, 1, 1, 0 },
	{  19,   0, EL_C, 2, 0, 0, 0 },
	{  38,  11, EL_C, 0, 0, 0, 1 },
	{  38, -11, EL_O, 0, 1, 1, 0 },
	{  57,  22, EL_O, 1, 1, 1, 0 }
};
static const nmr_peak_t pk_glycine[] = {
	{ 36, 2, M_S, 0, A_CH2 }, { 60, 3, M_S, 1, A_LABILES }
};

static const atom_t at_alanine[] = {
	{   0,  11, EL_C, 3, 1, 0, 0 },
	{  19,   0, EL_C, 1, 0, 0, 0 },
	{  19, -22, EL_N, 2, 1, 1, 0 },
	{  38,  11, EL_C, 0, 0, 0, 1 },
	{  38, -11, EL_O, 0, 1, 1, 0 },
	{  57,  22, EL_O, 1, 1, 1, 0 }
};
static const bond_t bd_alanine[] = {
	{ 0, 1, 1, 0, 0 }, { 1, 2, 1, 0, 0 }, { 1, 3, 1, 0, 0 },
	{ 3, 4, 2, 0, 0 }, { 3, 5, 1, 0, 0 }
};
static const nmr_peak_t pk_alanine[] = {
	{ 15, 3, M_D, 0, A_CH3 }, { 38, 1, M_Q, 0, A_NCH },
	{ 60, 3, M_S, 1, A_LABILES }
};

/* ================================================================== */
/* Table maitresse                                                     */
/*                                                                     */
/* Les isomeres d'une meme formule sont volontairement places cote a   */
/* cote : ils apparaissent dans cet ordre lors de la pagination.       */
/* ================================================================== */

#define MOL(a_, b_, c_, d_, nm, fam, at, bd, pk) \
	{ a_, b_, c_, d_, nm, fam, at, bd, pk, \
	  (uint8_t)(sizeof(at) / sizeof(atom_t)), \
	  (uint8_t)(sizeof(bd) / sizeof(bond_t)), \
	  (uint8_t)(sizeof(pk) / sizeof(nmr_peak_t)) }

/* Variante pour les molecules sans liaison a declarer (methane) */
static const bond_t bd_none[] = { { 0, 0, 1, 0, 0 } };

static const molecule_t DB[] = {
/* --- hydrocarbures --- */
MOL(1,  4, 0, 0, "Methane",        "Alcane",
    at_methane, bd_none, pk_methane),
MOL(2,  6, 0, 0, "Ethane",         "Alcane",
    at_ethane, bd_chain1, pk_ethane),
MOL(2,  4, 0, 0, "Ethene",         "Alcene",
    at_ethene, bd_double1, pk_ethene),
MOL(2,  2, 0, 0, "Ethyne",         "Alcyne",
    at_ethyne, bd_triple1, pk_ethyne),
MOL(6, 12, 0, 0, "Cyclohexane",    "Cycle sature",
    at_cyclohexane, bd_ring6, pk_cyclohexane),
MOL(6,  6, 0, 0, "Benzene",        "Aromatique",
    at_benzene, bd_benzene, pk_benzene),
MOL(7,  8, 0, 0, "Toluene",        "Aromatique",
    at_toluene, bd_sub1, pk_toluene),
MOL(8, 10, 0, 0, "Ethylbenzene",   "Aromatique",
    at_ethylbenzene, bd_sub2, pk_ethylbenzene),
MOL(8, 10, 0, 0, "para-Xylene",    "Aromatique para (2 signaux)",
    at_pxylene, bd_pxylene, pk_pxylene),
MOL(8, 10, 0, 0, "ortho-Xylene",   "Aromatique ortho",
    at_oxylene, bd_oxylene, pk_oxylene),
MOL(8,  8, 0, 0, "Styrene",        "Vinylaromatique",
    at_styrene, bd_styrene, pk_styrene),

/* --- alcools et ethers --- */
MOL(1,  4, 1, 0, "Methanol",       "Alcool primaire",
    at_methanol, bd_chain1, pk_methanol),
MOL(2,  6, 1, 0, "Ethanol",        "Alcool primaire",
    at_ethanol, bd_chain2, pk_ethanol),
MOL(2,  6, 1, 0, "Dimethylether",  "Ether-oxyde (1 seul signal)",
    at_dimethylether, bd_chain2, pk_dimethylether),
MOL(3,  8, 1, 0, "Propan-1-ol",    "Alcool primaire",
    at_propan1ol, bd_chain3, pk_propan1ol),
MOL(3,  8, 1, 0, "Propan-2-ol",    "Alcool secondaire",
    at_propan2ol, bd_iso3, pk_propan2ol),
MOL(4, 10, 1, 0, "Butan-1-ol",     "Alcool primaire",
    at_butan1ol, bd_chain4, pk_butan1ol),
MOL(4, 10, 1, 0, "tert-Butanol",   "Alcool tertiaire (2-methylpropan-2-ol)",
    at_tbutanol, bd_quat4, pk_tbutanol),
MOL(4, 10, 1, 0, "Ether diethylique", "Ether-oxyde",
    at_diethylether, bd_chain4, pk_diethylether),
MOL(6,  6, 1, 0, "Phenol",         "Phenol : OH tres variable",
    at_phenol, bd_sub1, pk_phenol),
MOL(7,  8, 1, 0, "Anisole",        "Ether aromatique",
    at_anisole, bd_sub2, pk_anisole),

/* --- aldehydes et cetones --- */
MOL(2,  4, 1, 0, "Ethanal",        "Aldehyde",
    at_ethanal, bd_carbonyl_end2, pk_ethanal),
MOL(3,  6, 1, 0, "Propanone",      "Cetone (1 seul signal)",
    at_propanone, bd_ketone3, pk_propanone),
MOL(3,  6, 1, 0, "Propanal",       "Aldehyde",
    at_propanal, bd_carbonyl_end3, pk_propanal),
MOL(4,  8, 1, 0, "Butanone",       "Cetone",
    at_butanone, bd_ketone4, pk_butanone),
MOL(4,  8, 1, 0, "Butanal",        "Aldehyde",
    at_butanal, bd_carbonyl_end4, pk_butanal),
MOL(7,  6, 1, 0, "Benzaldehyde",   "Aldehyde aromatique",
    at_benzaldehyde, bd_ar_carbonyl, pk_benzaldehyde),
MOL(8,  8, 1, 0, "Acetophenone",   "Cetone aromatique",
    at_acetophenone, bd_ar_ketone, pk_acetophenone),

/* --- acides et esters --- */
MOL(1,  2, 2, 0, "Acide methanoique", "Acide carboxylique",
    at_methanoic, bd_acid_head, pk_methanoic),
MOL(2,  4, 2, 0, "Acide ethanoique",  "Acide carboxylique",
    at_ethanoic, bd_acid2, pk_ethanoic),
MOL(2,  4, 2, 0, "Methanoate de methyle", "Ester (piege : isomere acide)",
    at_methylformate, bd_formate, pk_methylformate),
MOL(3,  6, 2, 0, "Acide propanoique", "Acide carboxylique",
    at_propanoic, bd_acid3, pk_propanoic),
MOL(3,  6, 2, 0, "Ethanoate de methyle", "Ester (2 singulets)",
    at_methylacetate, bd_ester4, pk_methylacetate),
MOL(4,  8, 2, 0, "Ethanoate d'ethyle", "Ester",
    at_ethylacetate, bd_ester5, pk_ethylacetate),
MOL(7,  6, 2, 0, "Acide benzoique",  "Acide aromatique",
    at_benzoic, bd_ar_acid, pk_benzoic),
MOL(8,  8, 2, 0, "Benzoate de methyle", "Ester aromatique",
    at_methylbenzoate, bd_ar_ester, pk_methylbenzoate),

/* --- composes azotes --- */
MOL(1,  5, 0, 1, "Methylamine",    "Amine primaire",
    at_methylamine, bd_chain1, pk_methylamine),
MOL(2,  7, 0, 1, "Ethylamine",     "Amine primaire",
    at_ethylamine, bd_chain2, pk_ethylamine),
MOL(6,  7, 0, 1, "Aniline",        "Amine aromatique",
    at_aniline, bd_sub1, pk_aniline),
MOL(2,  3, 0, 1, "Acetonitrile",   "Nitrile (DI = 2 pour C#N)",
    at_acetonitrile, bd_nitrile, pk_acetonitrile),
MOL(1,  3, 1, 1, "Formamide",      "Amide primaire",
    at_formamide, bd_acid_head, pk_formamide),
MOL(2,  5, 1, 1, "Acetamide",      "Amide primaire",
    at_acetamide, bd_acid2, pk_acetamide),
MOL(2,  5, 1, 1, "N-methylformamide", "Amide secondaire",
    at_nmethylformamide, bd_formate, pk_nmethylformamide),
MOL(6,  5, 2, 1, "Nitrobenzene",   "Derive nitre",
    at_nitrobenzene, bd_nitro, pk_nitrobenzene),
MOL(2,  5, 2, 1, "Glycine",        "Acide amine",
    at_glycine, bd_acid3, pk_glycine),
MOL(3,  7, 2, 1, "Alanine",        "Acide amine",
    at_alanine, bd_alanine, pk_alanine)
};

#define DB_SIZE ((int)(sizeof(DB) / sizeof(molecule_t)))

#endif /* CHEM_DB_H */
