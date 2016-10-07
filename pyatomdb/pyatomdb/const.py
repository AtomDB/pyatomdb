"""
This contains a list of constants, both physical and apec code related.

Version 0.1 - initial release
Adam Foster July 17th 2015
"""

#PHYSICAL CONSTANTS
KBOLTZ = 8.617385e-8#  /*!< in units of keV/K */
M_E = 2.7182818284590452354 # /*!< Euler e */
UPSILON_COLL_COEFF = 8.629e-6 # /*!< sqrt{2 pi / kB} hbar^2/m_e^{1.5} */
ERG_KEV = 1.60219e-9
ME_KEV = 510.0
RYDBERG = 0.013605804
BREMS_COEFF = 1.6542996e-20

HC_IN_KEV_A = 12.398425 #  hc for converting keV to And
HC_IN_ERG_A = 1.9862e-8  # /* = hc (erg A) = 12398.425 * 1.602e-12 */
SAF_COEFF = 3.30e-24

# APEC SPECIFIC CONSTANTS
MAX_UPS = 20
MAX_CHI = 200.0
MAX_IONREC = 20

E_UPSILON=1   # /*!< Electron upsilon values (unitless) */
E_RATE_COEFF=2 #/*!< Electron rate coefficient (cm^3/s) */
P_UPSILON=3    #/*!< Proton upsilon values (unitless) */
P_RATE_COEFF=4 #/*!< Proton rate coefficient (cm^3/s) */
EI_UPSILON=5   # /*!< Electron-impact ionization upsilon values (unitless) */

BURGESSTULLY=1 #/*!< Burgess-Tully-type data*/
CHIANTI_1=11   #/*!< CHIANTI pre-4.0 type 1 data (5 pt spline) */
CHIANTI_2=12   #/*!< CHIANTI pre-4.0 type 2 data (5 pt spline) */
CHIANTI_3=13   #/*!< CHIANTI pre-4.0 type 3 data (5 pt spline) */
CHIANTI_4=14   #/*!< CHIANTI pre-4.0 type 4 data (5 pt spline) */
CHIANTI_5=15   #/*!< CHIANTI pre-4.0 type 5 data (5 pt spline) */
CHIANTI_6=16   #/*!< CHIANTI pre-4.0 type 6 data (5 pt spline) */

CHIANTI4_1=21   #/*!< CHIANTI 4.0 type 1 data (9 pt spline) */
CHIANTI4_2=22   #/*!< CHIANTI 4.0 type 2 data (9 pt spline) */
CHIANTI4_3=23   #/*!< CHIANTI 4.0 type 3 data (9 pt spline) */
CHIANTI4_4=24   #/*!< CHIANTI 4.0 type 4 data (9 pt spline) */
CHIANTI4_5=25   #/*!< CHIANTI 4.0 type 5 data (9 pt spline) */
CHIANTI4_6=26   #/*!< CHIANTI 4.0 type 6 data (9 pt spline) */

SGC_1=31   #/*!< Sampson, Goett and Clark (1983) S-type He-like data */
SGC_2=32   #/*!< Sampson, Goett and Clark (1983) P-type He-like data */
SGC_3=33   #/*!< Sampson, Goett and Clark (1983) S-type H-like data */
KATO_NAKAZAKI_1=41  #/*!< Kato and Nakazaki (1989), ADNDT 42, 313 */
KATO_NAKAZAKI_2=42  #/*!< Kato and Nakazaki (1989), ADNDT 42, 313 */

#/* These must be spaced by at least MAX_UPS */

INTERP_E_UPSILON=100     #/*!< Include both left & right boundaries */
INTERP_P_UPSILON=200     #/*!< Include both left & right boundaries */
INTERP_E_RATE_COEFF=300  #/*!< Include both left & right boundaries */
INTERP_P_RATE_COEFF=400  #/*!< Include both left & right boundaries */

INTERP_E_UPS_OPEN=150    #/*!< Include neither boundary */
INTERP_P_UPS_OPEN=250    #/*!< Include neither boundary */
INTERP_E_RATE_OPEN=350   #/*!< Include neither boundary */
INTERP_P_RATE_OPEN=450   #/*!< Include neither boundary */

INTERP_E_UPS_INC_MIN=500   #/*!< Include only minimum; max is out */
INTERP_P_UPS_INC_MIN=600   #/*!< Include only minimum; max is out */
INTERP_E_RATE_INC_MIN=700  #/*!< Include only minimum; max is out */
INTERP_P_RATE_INC_MIN=800  #/*!< Include only minimum; max is out */

INTERP_E_UPS_INC_MAX=550   #/*!< Include only maximum; min is out */
INTERP_P_UPS_INC_MAX=650   #/*!< Include only maximum; min is out */
INTERP_E_RATE_INC_MAX=750  #/*!< Include only maximum; min is out */
INTERP_P_RATE_INC_MAX=850  #/*!< Include only maximum; min is out */

INTERP_I_UPSILON=900     #/*!< Include both left & right boundaries */

PROTON_BT=1001 #/*!< For Burgess-Tully Proton excitation rates */


CI_YOUNGER      =   51
EA_ARNROTH_LITHIUM =61
EA_ARNROTH         =62
EA_MAZZOTTA_IRON   =63
RR_SHULL  =71
RR_VERNER =72
RR_ARNRAY =73
RR_BADNELL =74
DR_MAZZOTTA =81
DR_BADNELL = 82
CI_DERE=100

INTERP_IONREC_RATE_COEFF =300
INTERP_IONREC_RATE_OPEN  =350
INTERP_IONREC_RATE_INC_MIN =700
INTERP_IONREC_RATE_INC_MAX =750


NO_PHOT= -1
HYDROGENIC= 0
CLARK= 1
VERNER= 2
XSTAR= 3
XSTAR_49_MIN =10000
XSTAR_49_MAX =19999
XSTAR_53_MIN =20000
XSTAR_53_MAX =29999


#Bremsstrahlung types:
HUMMER = 1
KELLOGG = 2
RELATIVISTIC = 3
BREMS_NONE = 4

#MISCfile types in FILEMAP files
ABUNDANCE = 10
HUMMER_TYPE = 11   
GAUNT_FF_TYPE = 13 


RRC_COEFF =  1.31e8
TOT_ABSACC = 1.e-21
TOT_RELACC = 0.0001
BIN_ABSACC = 1.e-21
BIN_RELACC = 0.01

TOLERANCE = 0.01

MIN_RRC_EXPONENT=-150.
MAX_RRC_EXPONENT=150.


ROMANIK   = 1
SAFRANOVA = 2

# ATOMDB ONLINE CONSTANTS
FILEDOWNLOAD = 1

# APEC SIZE LIMITS
NLEV_NOSPARSE = 5000
MIN_IONPOP = 1e-10
MIN_LEVPOP = 1e-40

