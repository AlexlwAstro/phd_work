This is a readme text file for imbeciles like me who are learning FRANEC/BaSTI for the first time!
This details the concepts underlined in the directory of the code.

***N.B.: FORTRAN 77 code is CASE INSENSITIVE!
***.2p1 files appear to be headers

CODE (_alw) folder - contains source code:
	atmosfer.2p1:
	Appears to be the header for atmospheric quantities. Atmospheric calculations are not considered after a **certain** threshold is reached by the model

	consts.2p1:
	Contains universal nature constants, such as Planck's constant, pi and G, together with astronomical unit conversions and shortcut labels for derived values, e.g. "mc2" = m*c*c (rest mass 		energy)

	maincom.2p1:
	Initialises parameters LIM (currently 18000), MAXPR (max number of***) and MAXNE (max number of elemental species to be considered) and makes COMMON various parameters used in different codes 	and routines

	diffus.f:
	FORTRAN source file containing diffusion-related subroutines:
		1. DIFFUS: calls in parameters from the equation of state to perform diffusion modelling.
		2. TBL: tabulates elemental abundances in diffusion via iteration loops.
		3. POLINO: carries out polynomial regression using matrices, which are also constructed here
		4. XMATR: carries out matrix multiplication of matrices X and Y and of X and X', its transpose.
		5. COEDIF: inverses the Burgers equations, written by Anne Thoul at Princeton.
		6. LUBKSB
		and
		7. LUDCMP are both used by COEDIF, and are from Numerical Recipes.

	energy.f:
	FORTRAN source file containing subroutines related to energy generation processes in the stellar interior:
		1. ENERGY: ****
		2. EPSI: **** N.B. Requires Centimetre-Gram-Second (CGS) units!!!
		3. EPSIG: ****
		4. NCROSS ****
		5. EVOLUT ****
		6. INDICIZZA ****
		7. RAP ****
		8. REANDET ****
		9. PLASMA ****
		10. SMALLRAP ****

	fisica.f:
	FORTRAN source file containing subroutines to calculate the equation of state (EOS) in detail:

	math.f:
	FORTRAN source file containing subroutines to calculate non-physical quantities, such as kernels and the ****

	mixing.f

	natmos.f

	optima.f

	main.f:
	FORTRAN source file combining everything into an evolution model**** via subroutines:
		*****

	

	modstart (initialise parameters here):

	1st field - model start point:
		1 = from modant
		2 = pre-main sequence
		3 = zero-age main sequence (ZAMS)
		4 = Zero-Age HB (ZAHB)
		5 = Horizontal Branch (HB)

	2nd field - number of layers in model (9999 = infinity)

	3rd field - number of steps between model saves (in case convergence fails, can retrace steps to last save)

	4th field - ?

	Binary fields following this - ?

FeHp006 (_alw) folder - contains data, both input and output, for the executable (.go) files:
	franec2p1p4.go:
	Executable file that runs the model calculations and stores the key results in grafi2p1

	bgt2p1.go:
	Executable file that tabulates the results from grafi2p1 into a more- (bigtab2p1) and less-detailed (asttab2p1) format

	If convection is turned on in modstart, the 

	flussi2p1.dat:
	Tabulates and stores the flux resulting from the individual fusion reaction types listed ****relative to the solar flux of these types.


