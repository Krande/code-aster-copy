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

subroutine rnomat(icesd, icesl, icesv, imap, nomcri, &
                  adrma, jtypma, k, optio, vala, &
                  valb, coefpa, nommat)
! person_in_charge: van-xuan.tran at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/rccome.h"
#include "asterfort/rcvale.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: icesd, icesl, icesv, imap, adrma, jtypma, k
    real(kind=8) :: vala, valb, coefpa
    character(len=8) :: nommat
    character(len=10) :: optio
    character(len=16) :: nomcri
! ----------------------------------------------------------------------
! BUT: RECUPERER POUR LA MAILLE COURANTE LE NOM DU MATERIAU DONNE PAR
!      L'UTILISATEUR.
! ----------------------------------------------------------------------
! ARGUMENTS :
!  ICESD    IN   I  : ADRESSE DE L'OBJET CHAM_ELEM_S.CESD.
!  ICESL    IN   I  : ADRESSE DE L'OBJET CHAM_ELEM_S.CESL.
!  ICESV    IN   I  : ADRESSE DE L'OBJET CHAM_ELEM_S.CESV.
!  IMAP     IN   I  : NUMERO DE LA MAILLE COURANTE.
!  NOMCRI   IN   K  : NOM DU CRITERE.
!  ADRMA*   IN   I  : ADRESSE DE LA MAILLE CORRESPONDANT AU NOEUD.
!  JTYPMA*  IN   I  : ADRESSE DU TYPE DE LA MAILLE.
!  K*       IN   I  : POINTEUR SERVANT AU TEST : MATERIAU UNIQUE OU NON.
!  OPTIO    IN   K  : CALCUL AUX POINTS DE GAUSS OU AUX NOEUDS.
!  VALA     OUT  R  : VALEUR DU COEFFICIENT A.
!  VALB     OUT  R  : VALEUR DU COEFFICIENT B.
!  COEFPA   OUT  R  : COEFFICIENT DE PASSAGE CISAILLEMENT-TRACTION.
!  NOMMAT   OUT  K  : NOM DU MATERIAU AFFECTE A LA MAILLE COURANTE.
!
! REMARQUE : * ==> VARIABLES N'AYANT DE SENS QUE DANS LE CAS DU CALCUL
!                  DE LA FATIGUE AUX NOEUDS.
! ----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: iad, numa, jtyp
    real(kind=8) :: r8b, v(1)
    integer(kind=8) :: icodre(1)
    character(len=8) :: ktyp, dimk, k8b
!     ------------------------------------------------------------------
!
!234567                                                              012
    r8b = 0.d0
    k8b = '        '
    call jemarq()
!
    if (optio .eq. 'DOMA_ELGA') then
!
! 1. RECUPERATION DU NOM DU MATERIAU AFFECTE A LA MAILLE COURANTE
!      _____________________________________________________
!     I                                                     I
!     I       CALCUL DE LA FATIGUE AUX POINTS DE GAUSS      I
!     I_____________________________________________________I
!
        call cesexi('C', icesd, icesl, imap, 1, &
                    1, 1, iad)
!C VERIFICATION HORS BORNES DEFINIES DANS CESMAT
!C OU COMPOSANTE NON AFFECTEE
        ASSERT(iad .gt. 0)
        nommat = zk8(icesv-1+iad)
!
    else if (optio .eq. 'DOMA_NOEUD') then
!
!      _____________________________________________________
!     I                                                     I
!     I           CALCUL DE LA FATIGUE AUX NOEUDS           I
!     I_____________________________________________________I
!
        numa = zi(adrma-1+imap)
        jtyp = jtypma-1+numa
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(jtyp)), ktyp)
        call dismoi('TYPE_TYPMAIL', ktyp, 'TYPE_MAILLE', repk=dimk)
!
        if (dimk .eq. 'VOLU') then
            k = k+1
            call cesexi('C', icesd, icesl, numa, 1, &
                        1, 1, iad)
!C VERIFICATION HORS BORNES DEFINIES DANS CESMAT
!C OU COMPOSANTE NON AFFECTEE
            ASSERT(iad .gt. 0)
            if (k .gt. 1) then
                if (nommat .ne. zk8(icesv-1+iad)) then
                    call utmess('F', 'FATIGUE1_33')
                end if
            else
                nommat = zk8(icesv-1+iad)
            end if
        end if
!
! CAS OU LES 1ERES MAILLES SCRUTEES NE SONT PAS VOLUMIQUE
! (ON PASSE A LA SUIVANTE)
        if (k .eq. 0) goto 999
!
    end if
!
! CAS CRITERE EST UNE FORMULE, ON NE RECUPERE QUE LE MATERIAU
    if (nomcri(1:7) .eq. 'FORMULE') then
        vala = 0.d0
        valb = 0.d0
        coefpa = 1.d0
        goto 999
    end if
!
    call rccome(nommat, 'CISA_PLAN_CRIT', icodre(1))
    if (icodre(1) .eq. 1) then
        call utmess('F', 'FATIGUE1_63')
    end if
!
! 2.1 RECUPERATION DES PARAMETRES ASSOCIES AU CRITERE MATAKE POUR
!     LA MAILLE COURANTE
!
    if (nomcri(1:14) .eq. 'MATAKE_MODI_AC') then
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'MATAKE_A', v(1), icodre(1), 0)
        vala = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_64')
        end if
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'MATAKE_B', v(1), icodre(1), 0)
        valb = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_65')
        end if
!
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'COEF_FLEX_TORS', v(1), icodre(1), 0)
        coefpa = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_66')
        end if
!
! 2.2 RECUPERATION DES PARAMETRES ASSOCIES AU CRITERE DE DANG VAN POUR
!     LA MAILLE COURANTE
!
    else if (nomcri(1:16) .eq. 'DANG_VAN_MODI_AC') then
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'D_VAN_A ', v(1), icodre(1), 0)
        vala = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_67')
        end if
!
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'D_VAN_B ', v(1), icodre(1), 0)
        valb = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_68')
        end if
!
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'COEF_CISA_TRAC', v(1), icodre(1), 0)
        coefpa = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_69')
        end if
    end if
!
! 2.3 RECUPERATION DES PARAMETRES ASSOCIES AU CRITERE MATAKE_MODI_AV
!     POUR LA MAILLE COURANTE
!
    if (nomcri(1:14) .eq. 'MATAKE_MODI_AV') then
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'MATAKE_A', v(1), icodre(1), 0)
        vala = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_70')
        end if
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'MATAKE_B', v(1), icodre(1), 0)
        valb = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_71')
        end if
!
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'COEF_FLEX_TORS', v(1), icodre(1), 0)
        coefpa = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_72')
        end if
    end if
!
! 2.4 RECUPERATION DES PARAMETRES ASSOCIES AU CRITERE DANG_VAN_MODI_AV
!     POUR LA MAILLE COURANTE
!
    if (nomcri(1:16) .eq. 'DANG_VAN_MODI_AV') then
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'D_VAN_A ', v(1), icodre(1), 0)
        vala = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_73')
        end if
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'D_VAN_B ', v(1), icodre(1), 0)
        valb = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_74')
        end if
!
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'COEF_CISA_TRAC', v(1), icodre(1), 0)
        coefpa = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_72')
        end if
    end if
!
! 2.5 RECUPERATION DES PARAMETRES ASSOCIES AU CRITERE FATEMI_SOCIE
!     POUR LA MAILLE COURANTE
!
    if (nomcri(1:16) .eq. 'FATESOCI_MODI_AV') then
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'FATSOC_A', v(1), icodre(1), 0)
        vala = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_75')
        end if
!
        valb = 1.0d0
!
        call rcvale(nommat, 'CISA_PLAN_CRIT', 0, k8b, [r8b], &
                    1, 'COEF_CISA_TRAC', v(1), icodre(1), 0)
        coefpa = v(1)
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_72')
        end if
!
    end if
!
999 continue
!
    call jedema()
!
end subroutine
