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

subroutine avcipr(nbvec, vectn, vectu, vectv, nbordr, &
                  kwork, sommw, vwork, tdisp, tspaq, &
                  ipgn, nomcri, nomfor, fordef, fatsoc, &
                  proaxe, pseuil, method, ncycl, jvmin, &
                  jvmax, jomin, jomax)
! aslint: disable=W1306,W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
!
#include "asterfort/avenca.h"
#include "asterfort/avpeak.h"
#include "asterfort/avpic2.h"
#include "asterfort/avrain.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/projax.h"
#include "asterfort/propla.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbvec, nbordr, kwork
    integer(kind=8) :: sommw, tdisp, tspaq, ipgn
    real(kind=8) :: vectn(3*nbvec), vectu(3*nbvec), vectv(3*nbvec)
    real(kind=8) :: vwork(tdisp), fatsoc, pseuil
    character(len=16) :: nomcri, nomfor, proaxe
    character(len=8) :: method
    integer(kind=8) :: ncycl(nbvec)
    integer(kind=8) :: jomin, jomax, jvmin, jvmax
!    integer :: omin(nbvec*(nbordr+2)), omax(nbvec*(nbordr+2))
!    real(kind=8) :: vmin(nbvec*(nbordr+2)), vmax(nbvec*(nbordr+2))
    aster_logical :: fordef
!
! BUT:  POUR LA FATIGUE A AMPLITUDE VARIABLE
!       A PARTIR DE l'HISTOIRE DE CISSAILLEMENT, PROJETER SUR UN OU
!       2 AXES ET DETERMINER DES PICS PAR LE COMPTAGE DE RAINFLOW
! REMARQUE: CETTE SUBROUTINE EST APPLICABLE POUR UN NOEUD OU UN POINT
!           GAUSSE
!
! ----------------------------------------------------------------------
! ARGUMENTS :
!  NBVEC    IN  I  : NOMBRE DE VECTEURS NORMAUX.
!  VECTN    IN  R  : VECTEUR CONTENANT LES COMPOSANTES DES
!                    VECTEURS NORMAUX.
!  VECTU    IN  R  : VECTEUR CONTENANT LES COMPOSANTES DES
!                    VECTEURS u DU PLAN DE CISAILLEMENT.
!  VECTV    IN  R  : VECTEUR CONTENANT LES COMPOSANTES DES
!                    VECTEURS v DU PLAN DE CISAILLEMENT.
!  NBORDR   IN  I  : NOMBRE DE NUMEROS D'ORDRE.
!  KWORK    IN  I  : KWORK = 0 ON TRAITE LA 1ERE MAILLE DU PAQUET
!                              MAILLES OU LE 1ER NOEUD DU PAQUET DE
!                              NOEUDS;
!                    KWORK = 1 ON TRAITE LA IEME (I>1) MAILLE DU PAQUET
!                              MAILLES OU LE IEME NOEUD DU PAQUET
!                              DE NOEUDS.
!  SOMMW    IN  I  : SOMME DES POINTS DE GAUSS OU DES NOEUDS DES N
!                    MAILLES PRECEDANT LA MAILLE COURANTE.
!  VWORK    IN  R  : VECTEUR DE TRAVAIL CONTENANT
!                    L'HISTORIQUE DES TENSEURS DES CONTRAINTES
!                    ATTACHES A CHAQUE POINT DE GAUSS OU NOEUD DES
!                    MAILLE OU NOEUD DU <<PAQUET>> DE MAILLES OU
!                    DE NOEUDS.
!  TDISP    IN  I  : DIMENSION DU VECTEUR VWORK
!  TSPAQ    IN  I  : TAILLE DU SOUS-PAQUET DU <<PAQUET>> DE MAILLES
!                    OU DE NOEUDS COURANT.
!  I        IN  I  : IEME POINT DE GAUSS OU IEME NOEUD.
!  NOMCRI   IN  K16: NOM DU CRITERE D'ENDOMMAGEMENT PAR FATIGUE.
!  FATSOC   IN  R  : COEFFICIENT PERMETTANT D'UTILISER LES MEMES
!                    ROUTINES POUR LE TRAITEMENT DES CONTRAINTES ET
!                    DES DEFORMATIONS.
!  PROAXE    IN   K16: TYPE DE PROJECTION (UN OU DEUX AXES).
!  PSEUIL    IN   R  : SEUIL DE DTECTION DES PICS
!
!  METHOD    IN   K  : METHODE D'EXTRACTION DES PICS, PAR EXEMPLE :
!                     RAINFLOW.
! NCYCL     OUT  I  : NOMBRE DE CYCLES ELEMENTAIRES POUR TOUS LES
!                     VECTEURS NORMAUX.
! JVMIN      OUT  I  : ADDRESEE JEUVEUX DES VALEURS MIN DES CYCLES ELEMENTAIRES
!                     POUR TOUS LES VECTEURS NORMAUX.
! JVMAX      OUT  I  : ADDRESEE JEUVEUX DES VALEURS MAX DES CYCLES ELEMENTAIRES
!                     POUR TOUS LES VECTEURS NORMAUX.
! JOMIN      OUT  I  : ADDRESEE JEUVEUX DES NUMEROS D'ORDRE ASSOCIES AUX
!                      VALEURS MIN DESCYCLES ELEMENTAIRES POUR TOUS LES
!                      VECTEURS  NORMAUX.
! JOMAX      OUT  I  : ADDRESEE JEUVEUX DES NUMEROS D'ORDRE ASSOCIES AUX
!                      VALEURS MAX DES CYCLES ELEMENTAIRES POUR TOUS LES
!                      VECTEURS NORMAUX
! REMARQUE : CETTE ROUTINE SERT POUR LE TRAITEMENT DES POINTS DE GAUSS
!            ET DES NOEUDS.
! ----------------------------------------------------------------------
!
! VECTRA  VECTEUR DE TRAVAIL CONTENANT
!         LES COMPOSANTES u ET v DU VECTEUR TAU
!         (CONTRAINTE DE CISAILLEMENT) OU
!         GAMMA (DEFORMATION DE CISAILLEMENT), POUR TOUS LES
!         NUMEROS D'ORDRE DE CHAQUE VECTEUR NORMAL.
! LSIG0   VARIABLE LOGIQUE QUI INDIQUE :
!         - LSIG0 = FALSE --> CAS GENERAL, LES CONTRAINTES
!                             SONT DIFFERENTES DE ZERO ;
!         - LSIG0 =  TRUE --> LES CONTRAINTES SONT NULLES
!                             A TOUS LES PAS DE TEMPS, QUEL
!                             QUE SOIT LE VECTEUR NORMAL.
! IFLAG  VECTEUR DE DRAPEAUX QUI INDIQUE :
!         - IFLAG(i) = 0 --> CAS GENERAL ;
!         - IFLAG(i) = 1 --> CAS OU LES POINTS DANS LE
!                            PLAN DE CISAILLEMENT SONT
!                            ALIGNES VERTICALEMENT ;
!         - IFLAG(i) = 2 --> CAS OU LES POINTS DANS LE
!                            PLAN DE CISAILLEMENT SONT
!                            ALIGNES HORIZONTALEMENT ;
!         - IFLAG(i) = 3 --> CAS OU LES POINTS DANS LE
!                            PLAN DE CISAILLEMENT SONT
!                            CONTENUS DANS UN CADRE DE
!                            COTES INFERIEURS A EPSILO.
! RMIMA  VECTEUR CONTENANT LES COORDONNEES DES POINTS
!        EXTREMES DU CADRE (CUMIN, CUMAX, CVMIN, CVMAX)
!        POUR TOUS LES VECTEURS NORMAUX.
! RAXE   VECTEUR CONTENANT L'AMPLITUDE DES POINTS
!        PROJETES.
! NPOIN  NOMBRE DE PICS DETECTES POUR TOUS LES VECTEURS
!        NORMAUX.
! VALPOI VALEUR DES PICS DETECTES POUR TOUS LES VECTEURS
!        NORMAUX.
! VALORD NUMEROS D'ORDRE ASSOCIES AUX PICS DETECTES POUR
!        TOUS LES VECTEURS NORMAUX.
! NPIC   NOMBRE DE PICS DETECTES POUR TOUS LES VECTEURS
!        NORMAUX APRES REARANGEMENT DES PICS.
! PIC    VALEUR DES PICS DETECTES POUR TOUS LES VECTEURS
!        NORMAUX APRES REARANGEMENT DES PICS.
! ORDPIC NUMEROS D'ORDRE ASSOCIES AUX PICS DETECTES POUR
!        TOUS LES VECTEURS NORMAUX APRES REARANGEMENT
!        DES PICS.
! RTRV   VECTEUR DE TRAVAIL REEL (POUR LES POINTS)
! ITRV   VECTEUR DE TRAVAIL ENTIER (POUR LES NUME_ORDRE)
!
!    real(kind=8) :: vectra(2*nbvec*nbordr), rmima(4*nbvec)
!    integer :: iflag(nbvec), itrv(2*(nbordr+2))
!    aster_logical :: lsig0
!    real(kind=8) :: raxe(nbvec*nbordr), valpoi(nbvec*nbordr)
!    integer :: npoin(nbvec), valord(nbvec*nbordr)
!    integer :: npic(nbvec), ordpic(nbvec*(nbordr+2))
!    real(kind=8) :: pic(nbvec*(nbordr+2)), rtrv(nbordr+2)
!
    real(kind=8) :: rmima(4*nbvec)
    integer(kind=8) :: iflag(nbvec)
    aster_logical :: lsig0
    integer(kind=8) :: npoin(nbvec)
    integer(kind=8) :: npic(nbvec)
!
    integer(kind=8) :: jvectr, jitrv, jraxe, jvalpo, jvalor, jordpi, jpic, jrtrv
!
!      REAL*8        CUDOMX, NXM, NYM, NZM
!     ------------------------------------------------------------------
!
!  PROJECTION DE L'HISTORIQUE DU CISAILLEMENT DANS UN PLAN
    call jemarq()
!
! Real
    call wkvect('&&AVCIPR_VECTRA', 'V V R', 2*nbvec*nbordr, jvectr)
    call wkvect('&&AVCIPR_RAXE', 'V V R', nbvec*nbordr, jraxe)
    call wkvect('&&AVCIPR_VALPOI', 'V V R', nbvec*nbordr, jvalpo)
    call wkvect('&&AVCIPR_PIC', 'V V R', nbvec*(nbordr+2), jpic)
    call wkvect('&&AVCIPR_RTRV', 'V V R', (nbordr+2), jrtrv)
!
! Integer
    call wkvect('&&AVCIPR_ITVR', 'V V I', 2*(nbordr+2), jitrv)
    call wkvect('&&AVCIPR_VALORD', 'V V I', nbvec*nbordr, jvalor)
    call wkvect('&&AVCIPR_ORPIC', 'V V I', nbvec*(nbordr+2), jordpi)
!
    call propla(nbvec, vectn, vectu, vectv, nbordr, &
                kwork, sommw, vwork, tdisp, tspaq, &
                ipgn, nomcri, nomfor, fordef, fatsoc, &
                jvectr)
!
! CALCUL DU DOMMAGE MAX ET DU VECTEUR NORMAL ASSOCIE POUR
! LE NOEUD/POINT GAUSS COURANT DE LA MAILLE COURANTE.
!
! 1. REMISE A ZERO DU VECTEUR DE TRAVAIL CONTENANT LES VALEURS DE
!    DELTA_TAU POUR UN NOEUD ET DU VECTEUR DE TRAVAIL
!    PERMETTANT DE POINTER SUR LE VECTEUR NORMAL ASSOCIE.
!
! 2. ENCADREMENT DES POINTS DANS LE PLAN
!
    lsig0 = .false.
!
    call avenca(jvectr, nbvec, nbordr, lsig0, iflag, &
                rmima)
!
!       IF (LSIG0) THEN
!          CUDOMX = 0.0D0
!          NXM = 0.0D0
!          NYM = 0.0D0
!          NZM = 1.0D0
!          GOTO 555
!       ENDIF
!
! 3. PROJECTION DE L'HISTORIQUE DE CHARGEMENT SUR UN OU DEUX AXES
!
    call projax(jvectr, nbvec, nbordr, proaxe, iflag, &
                rmima, jraxe)
!
! 4. COMPTAGE RAINFLOW (NORME AFNOR + POSTDAM)
!
! 4.1 PREMIER FILTRAGE DES PICS DE LA FONCTION
!
    call avpeak(jraxe, nbvec, nbordr, pseuil, iflag, &
                npoin, jvalpo, jvalor)
!
! 4.2 REARANGEMENT ET EXTRACTION DES PICS
!
!
    call avpic2(method, nbvec, nbordr, jrtrv, jitrv, &
                npoin, jvalpo, jvalor, npic, jpic, &
                jordpi)
!
! 4.3 COMPTAGE RAINFLOW
!
!
    call avrain(nbvec, nbordr, jitrv, npic, jpic, &
                jordpi, fatsoc, ncycl, jvmin, jvmax, &
                jomin, jomax)
!
!
    call jedetr('&&AVCIPR_VECTRA')
    call jedetr('&&AVCIPR_ITVR')
    call jedetr('&&AVCIPR_RAXE')
    call jedetr('&&AVCIPR_VALPOI')
    call jedetr('&&AVCIPR_VALORD')
    call jedetr('&&AVCIPR_ORPIC')
    call jedetr('&&AVCIPR_PIC')
    call jedetr('&&AVCIPR_RTRV')
!
    call jedema()
!
end subroutine
