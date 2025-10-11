! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
! 
! -------------------------------------------------------------------------
!

! constantes
! constante des gaz parfait en J/K/mol
#define KGAZP 8.314d0
!   énergie d'activation de la viscosté de l'eau en J/K/mol
#define EASURRW 2059.d0
! zero Kelvin
#define ZEROKLV 273.15d0
! 80 degré
#define T80DEG 80.d0
! tolérance
#define EPSIL 1.d-6
! volume molaire de l’ettringite différé
#define VMAFT 715.d-6
#define VMAFM 254.6d-6
! Multiplicateur non linéaire maximum de potentiel de fluage
#define XMAX 25.d0
! water molar mass (g/mol)
#define WATERMOLARMASS 18.01528d-3
! atm to Pa
#define AtmToPa 101325.d0

! pour indices dans xmat

#define YOUN  1
#define NU    2  

#define HYDR  5  
#define HYDS  6  
#define RT    7  
#define REF   8  
#define RC    9 
#define DELT  10
#define BETA  11
#define EPT   12
#define HRGI  13
#define VVRG  14
#define KGEL  15
#define GFT   16
#define EKFL  17
#define YKSY  18
#define XFLU  19
#define TAUK  20
#define TAUM  21
#define NRJM  22
#define DT80  23
#define BSHR  24
#define MSHR  25
#define PORO  26
#define VRAG  27
#define NRJF  28
#define SFLD  29
#define MVGN  30
#define EPC   31
#define EKDC  32
#define EKRG  33
#define GFR   34
#define ALAT  35
#define KRGI  36
#define TREF  37
#define TSTH  38
#define DFMX  39
#define TAUG  40
#define NRJG  41
#define SRSG  42
#define TRAG  43
#define DIM3  44
#define TDEF  45
#define NRJP  46
#define SRSD  47
#define VDEF  48
#define CNAD  49
#define SSAD  50
#define CNAK  51
#define CNAB  52
#define EXND  53
#define EXMD  54
#define TTDD  55
#define TDID  56
#define TFID  57
#define NRJD  58
#define TTRD  59
#define TTKF  60
#define HPEV  61
           
#define YOUM  62
#define NUM   63
#define NREN  64




! nombre de variables internes concernant le beton uniquement
#define NVARBE 114
! nombre de variables internes par renfort
#define NVARENF 7

! pour indices dans var0 et varf (variables internes)

! deformations elastiques
#define EPE(i) i
! deformation de l etage de Kelvin
#define EPK(i) 7-1+i 
! deformations de l etage de Maxwell          
#define EPM(i) 13-1+i
! contraintes effective squelette solide (sans pression rgi)          
#define SIG(i) 19-1+i
! dissipation etage de Mawell          
#define PHIM 25
! energie elastique          
#define WELA 26
! endomagement par fluage
#define DFLU 27
! variation de volume meca total
#define TEPS 28
! variation de volume par fissuration rgi          
#define TEPG 29
! deformation plastique de traction         
#define EPTi(i) 30-1+i
! de36formation plastique fissuration rgi
#define EPGi(i) 36-1+i
! deformations plastiques de druker prager       
#define EPCi(i) 42-1+i
! hydratation fin de pas
#define HYDF 48
! endommagement thermique isotrope
#define DTHE 49
! contraintes elastiques de l etage de Kelvin
#define SKE(i) 50-1+i
! pression capillaire de l eau
#define PSHR 56
! perte de viscosité par séchage (coeff de consolidation par séchage)
#define CSHR 57
! volume d eau capillaire (necessaire sous iteration fluage3d)
#define WSHR 58
! vide pour l eau  capillaire (necessaire sous iteration fluage3d)
#define VSHR 59
! potentiel des phases neoformees
#define PHIG 60
! pression de RGI
#define PRGI 61
! avancement AAR
#define AAAR 62
! avancement DEF
#define ADEF 63
! indicateur premier pas (1 si premier pas passé sinon 0)
#define PPAS 64
! coeff de Biot effectif pour les RGI
#define BIOG 65
! coeff de Biot effectif pour la capillarite
#define BIOW 66
! deformation plastique equivalente en cisaillement
#define EPLC 67
! deformation plastique maximale atteinte en traction
#define EMT(i) 68-1+i
! Cumul des surpessions capillaires en dessiccation sous charge
#define DSW(i) 74-1+i
! endommagements principaux de traction directe 
#define DTL(i) 80-1+i
! endommagements principaux de refermeture 
#define DRL(i) 83-1+i
! endommagements principaux de traction diffus RGI
#define DTG(i) 86-1+i
! endommagements principaux de compression difffus RGI
#define DCG(i) 89-1+i
! partie plastique de l ouverture de fissure localisee
#define WL(i) 92-1+i
! endommagement de compression
#define DC 95
! endommagement de traction diffus pre pic
#define DTPP 96
! sulfo aluminates pour la def 
#define AFT1 97
#define AFM1 98
#define AFT2 99
#define AFM2 100
#define ATIL 101
#define STIL 102
! tenseur des ouvertures plastiques de fissures en base globale
#define WID(i) 103-1+i
! erreur commise sur la dissipation d'énergie de fissuration en traction
#define ERGF 109
! endommagement de traction global 
#define DT0 110
! Ouverture de fissure maximale 
#define WPL0 111
! Endommagement de traction de RGI global 
#define DTG0 112
! Endommagement de compression de RGI global 
#define DCG0 113
! Pression de gel maximal atteinte 
#define PGMAX 114
! variables pour les renforts acier (5 directions possibles)
#define EPSIN(i) NVARBE+NVARENF*(i-1)+1
#define EPSEQ(i) NVARBE+NVARENF*(i-1)+2
#define SNR(i) NVARBE+NVARENF*(i-1)+3
#define ERM(i) NVARBE+NVARENF*(i-1)+4
#define MUR(i) NVARBE+NVARENF*(i-1)+5
#define ERK(i) NVARBE+NVARENF*(i-1)+6
#define SPR(i) NVARBE+NVARENF*(i-1)+7
! Contrainte dans la matrice seule
#define SBE(i) 150-1+i
! Contrainte principale de la matrice de beton 
#define SPM1 156
#define SPM2 157
#define SPM3 158
