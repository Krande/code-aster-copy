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

subroutine mpjeft(corres)
!
!   MPJEFT : CREATION TABLE DE CORRESPONDANCE POUR PROJ_MESU_MODAL
!
!   IN/OUT : CORRES : NOM DE LA TABLE
!
! ----------------------------------------------------------------------
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/pacoa2.h"
#include "asterfort/pj2dco.h"
#include "asterfort/pj3dco.h"
#include "asterfort/pj4dco.h"
#include "asterfort/pj5dco.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    character(len=16) :: corres
!
!
!
!
    character(len=4) :: cdim1, exivol, exipou, exirdm, exipla, exicoq
    character(len=8) :: noma1, noma2, model1, model2, labk8
    character(len=8) :: lisin1, lisin2, lisou1, lisou2
    character(len=16) :: tymocl(2), motcle(2)
    integer(kind=8) :: ndim, ncas, n1, nbocc, iocc, nbno2, nuno1, nuno2
    integer(kind=8) ::  idecal, ino, irep
    integer(kind=8) :: llin1, llin2, llou2, nbncal, nbnlis
    integer(kind=8) :: kk, nbnmes, nbno1, jxxk1, iaconb, iaconu, iacocf, i
    integer(kind=8) :: nbold
    integer(kind=8) :: ndecal, it1, it2, it3, ifres
    real(kind=8) :: coef, rbid
    integer(kind=8), pointer :: linonu1(:) => null()
    integer(kind=8), pointer :: linonu2(:) => null()
    aster_logical :: l_dmax
!----------------------------------------------------------------------
!       DESCRIPTION DE LA SD CORRESP_2_MAILLA : CORRES
!       --------------------------------------------------------
!          (CORRESPONDANCE ENTRE LES 2 MODELES MODEL1 ET MODEL2)
!  DESCRIPTION DE LA SD CORRESP_2_MAILLA :
!
!  CORRESP_2_MAILLA (K16)  ::= RECORD
!     '.PJXX_K1'  : S V K24 LONG=5
!     '.PJEF_NB'  : S V I  LONG=NNO2 (= NB_NO(M2))
!     '.PJEF_NU'  : S V I  LONG=LONT
!     '.PJEF_CF'  : S V R  LONG=LONT
!
!
!     '.PJXX_K1' (1) : NOM DU MAILLAGE 1 : M1
!     '.PJXX_K1' (2) : NOM DU MAILLAGE 2 : M2
!
!     '.PJEF_NB' (INO2) : NOMBRE DE NOEUDS DE M1 QUI DOIVENT SERVIR
!                          A L'INTERPOLATION DU NOEUD INO2 DE M2
!
!     '.PJEF_NU' : CONTIENT LES NUMEROS DES NOEUDS DE M1 SERVANT A
!                  L'INTERPOLATION DES NOEUDS DE M2 (MIS BOUT A BOUT)
!     '.PJEF_CF' : CONTIENT LES COEFFICIENTS POUR LES NOEUDS DE
!                  M1 SERVANT A L'INTERPOLATION DES NOEUDS DE M2
!                  (MIS BOUT A BOUT)
!
!    EXEMPLE D'UTILISATION :
!      ON VEUT SAVOIR COMMENT INTERPOLER INO2 A PARTIR DU MAILLAGE M1
!      SOIT NBNO1='.PJEF_NB'(INO2)
!      SOIT DECAL= SOMME POUR INO<INO2 DE '.PJEF_NB'(INO)
!      VAL2(INO2)=0
!      DO I=1,NBNO1
!        NUNO1='.PJEF_NU' (DECAL+I)
!        COEFR='.PJEF_CF' (DECAL+I)
!        VAL2(INO2)=VAL2(INO2)+COEFR*VAL1(UNO1)
!      END DO
!----------------------------------------------------------------------
! DEB ------------------------------------------------------------------
    call jemarq()
!
    call getvid('MODELE_CALCUL', 'MODELE', iocc=1, scal=model1, nbret=n1)
    call getvid('MODELE_MESURE', 'MODELE', iocc=1, scal=model2, nbret=n1)
!
    call dismoi('NOM_MAILLA', model1, 'MODELE', repk=noma1)
    call dismoi('NOM_MAILLA', model2, 'MODELE', repk=noma2)
!
    call dismoi('NB_NO_MAILLA', noma2, 'MAILLAGE', repi=nbnmes)
!
    call dismoi('NB_NO_MAILLA', noma1, 'MAILLAGE', repi=nbncal)
!
!     DETERMINATION DE LA DIMENSION DE L'ESPACE (NDIM) :
!     --------------------------------------------------------
    call dismoi('Z_CST', noma1, 'MAILLAGE', repk=cdim1)
!
    if (cdim1 .eq. 'OUI') then
        ndim = 2
    else
        ndim = 3
    end if
!
!     PAR DEFAUT ON FAIT UNE INTERPOLATION SUR TOUS LES NOEUDS
!     DU MAILLAGE 2 (MAILLAGE EXPERIMENTAL)
!     ON SURCHARGE ENSUITE LES CORRESPONDANCES MANUELLES
!
!     DETERMINATION DU CAS DE FIGURE : 2D, 3D OU 2.5D : NCAS
!     --------------------------------------------------------
    if (ndim .eq. 2) then
        ncas = 2
    else if (ndim .eq. 3) then
        call dismoi('EXI_ELTVOL', model1, 'MODELE', repk=exivol)
        if (exivol .eq. 'OUI') then
            ncas = 3
        else
            ncas = 4
        end if
    end if
    call dismoi('EXI_RDM', model1, 'MODELE', repk=exirdm)
    irep = 0
    call dismoi('EXI_PLAQUE', model1, 'MODELE', repk=exipla)
    call dismoi('EXI_COQ3D', model1, 'MODELE', repk=exicoq)
    if (exipla .eq. 'OUI') irep = 1
    if (exicoq .eq. 'OUI') irep = 1
    if ((exirdm .eq. 'OUI') .and. (irep .eq. 0)) ncas = 0
    call dismoi('EXI_POUX', model1, 'MODELE', repk=exipou)
    if (exipou .eq. 'OUI') ncas = 5
!
!
! POUR L'INSTANT ON SE LIMITE AU CAS OU IL N'Y A PAS DE MELANGE
!  ELEMENT VOLUMIQUE / ELEMENT SURFACIQUE / ELEMENT LINEIQUE
!
    l_dmax = ASTER_FALSE
    if (ncas .eq. 2) then
        call pj2dco('TOUT', model1, model2, 0, [0], &
                    0, [0], ' ', ' ', corres, &
                    l_dmax, rbid, 0.d0)
    else if (ncas .eq. 3) then
        call pj3dco('TOUT', model1, model2, 0, [0], &
                    0, [0], ' ', ' ', corres, &
                    l_dmax, rbid, 0.d0)
    else if (ncas .eq. 4) then
        call pj4dco('TOUT', model1, model2, 0, [0], &
                    0, [0], ' ', ' ', corres, &
                    l_dmax, rbid, 0.d0)
    else if (ncas .eq. 5) then
        call pj5dco(model1, model2, corres)
    else
!
        call utmess('A', 'ALGORITH6_19')
!
!     CREATION DE LA SD CORRESP_2_MAILLA : CORRES (PROCHE)
!     ---------------------------------------------------
        call wkvect(corres//'.PJXX_K1', 'V V K24', 5, jxxk1)
        call wkvect(corres//'.PJEF_NB', 'V V I', nbnmes, iaconb)
        call wkvect(corres//'.PJEF_NU', 'V V I', nbnmes, iaconu)
        call wkvect(corres//'.PJEF_CF', 'V V R', nbnmes, iacocf)
!
        zk24(jxxk1-1+1) = noma1
        zk24(jxxk1-1+2) = noma2
        zk24(jxxk1-1+3) = 'ELEM'
!
        do ino = 1, nbnmes
            zi(iaconb-1+ino) = 1
            zr(iacocf-1+ino) = 1.d0
        end do
!
! CREATION DES LISTES DE NOEUDS UTILISEES
! ***************************************
!
!       -> LISTE DES NOEUDS DE MESURE
        lisin1 = 'LINNOMES'
!
!       -> LISTE DES NOEUDS DU MODELE NUMERIQUE
        lisin2 = 'LINNOCAL'
!
!       -> LISTE TRIEE DES NOEUDS DE MESURE
        lisou1 = 'LOUNOMES'
!
!       -> LISTE TRIEE DES NOEUDS NUMERIQUES EN VIS-A-VIS DE LISOU1
        lisou2 = 'LOUNOCAL'
!
! ALLOCATION ET REMPLISSAGE DES LISTES DE NOEUDS DE MESURE
! ********************************************************
        call wkvect(lisin1, 'V V K8', nbnmes, llin1)
!
        do ino = 1, nbnmes
            zk8(llin1-1+ino) = int_to_char8(ino)
        end do
!
! ALLOCATION ET REMPLISSAGE DE LA LISTE DES NOEUDS NUMERIQUES
! ***********************************************************
        call wkvect(lisin2, 'V V K8', nbncal, llin2)
!
        do ino = 1, nbncal
            zk8(llin2-1+ino) = int_to_char8(ino)
        end do
!
! RECHERCHE DES NOEUDS EN VIS-A-VIS
! *********************************
        call pacoa2(lisin1, lisin2, nbnmes, nbncal, noma2, &
                    noma1, lisou1, lisou2, nbnlis)
!
        call jeveuo(lisou2, 'L', llou2)
!
        if (nbnlis .ne. nbnmes) then
            call utmess('F', 'ALGORITH6_20')
        end if
!
        do ino = 1, nbnmes
            zi(iaconu-1+ino) = char8_to_int(zk8(llou2-1+ino))
        end do
!
        call jedetr(lisin1)
        call jedetr(lisin2)
        call jedetr(lisou1)
        call jedetr(lisou2)
!
    end if
!
!     SURCHARGE SI CORRESPONDANCE MANUELLE
!
    call getfac('CORR_MANU', nbocc)
    if (nbocc .gt. 0) then
!
!       CORRESPONDANCE MANUELLE SUR CERTAINS NOEUDS
!       ------------------------
!
        do iocc = 1, nbocc
!
!        -- RECUPERATION DE LA LISTE DE NOEUDS LNO1 (CALCUL):
!        ----------------------------------------------
            motcle(1) = 'NOEU_CALCUL'
            tymocl(1) = 'NOEUD'
            call reliem(' ', noma1, 'NU_NOEUD', 'CORR_MANU', iocc, &
                        1, motcle, tymocl, '&&PJEFTE.LINONU1', nbno1)
            if (nbno1 .eq. 0) then
                call utmess('F', 'ALGORITH6_21')
            end if
            call jeveuo('&&PJEFTE.LINONU1', 'L', vi=linonu1)
!
!
!        -- RECUPERATION DE LA LISTE DE NOEUDS LNO2 (MESURE):
!        ----------------------------------------------
            motcle(1) = 'NOEU_MESURE'
            tymocl(1) = 'NOEUD'
            call reliem(' ', noma2, 'NU_NOEUD', 'CORR_MANU', iocc, &
                        1, motcle, tymocl, '&&PJEFTE.LINONU2', nbno2)
            if (nbno2 .eq. 0) then
                call utmess('F', 'ALGORITH6_22')
            end if
            call jeveuo('&&PJEFTE.LINONU2', 'L', vi=linonu2)
!
!        -- REACTUALISATION DU CORRESP_2_MAILLA POUR IOCC
!        ----------------------------------------------
            if (nbno1 .ne. nbno2) then
                call utmess('F', 'ALGORITH6_23')
            end if
            if (nbno1 .gt. 1) then
                call utmess('F', 'ALGORITH6_24')
            end if
!        -- RECUPERATION DES NUMEROS DES NOEUDS
            nuno1 = linonu1(1)
            nuno2 = linonu2(1)
!
!        -- RECUPERATION DES ADRESSES SD CORRES
!     ------------------------------------------------
            call jeveuo(corres//'.PJEF_NB', 'E', iaconb)
            call jeveuo(corres//'.PJEF_NU', 'E', iaconu)
            call jeveuo(corres//'.PJEF_CF', 'E', iacocf)
!
! RECUPERATION DE LA DIMENSION DE PJEF_NU OU PJEF_CF
            idecal = 0
            do i = 1, nbnmes
                idecal = idecal+zi(iaconb-1+i)
            end do
! CREATION DES VECTEURS TAMPON
            call wkvect('TAMPON1', 'V V I', nbnmes, it1)
            call wkvect('TAMPON2', 'V V I', idecal, it2)
            call wkvect('TAMPON3', 'V V R', idecal, it3)
! RECOPIE DE PJEF_NB PJEF_NU ET PJEF_CF
            do i = 1, nbnmes
                zi(it1-1+i) = zi(iaconb-1+i)
            end do
            do i = 1, idecal
                zi(it2-1+i) = zi(iaconu-1+i)
                zr(it3-1+i) = zr(iacocf-1+i)
            end do
! DESTRUCTION DES ANCIENS PJEF_NB PJEF_NU ET PJEF_CF
            call jedetr(corres//'.PJEF_NB')
            call jedetr(corres//'.PJEF_NU')
            call jedetr(corres//'.PJEF_CF')
!
!        -- DECALAGE DES DONNEES
            idecal = 0
            do i = 1, nuno2-1
                idecal = idecal+zi(it1-1+i)
            end do
!
            nbold = zi(it1-1+nuno2)
            zi(it1-1+nuno2) = 1
            zi(it2-1+idecal+1) = nuno1
            zr(it3-1+idecal+1) = 1.d0
!
            ndecal = 0
            do i = nuno2+1, nbnmes
                ndecal = ndecal+zi(it1-1+i)
            end do
!
            do i = 1, ndecal
                zi(it2-1+idecal+1+i) = zi(it2-1+idecal+nbold+i)
                zr(it3-1+idecal+1+i) = zr(it3-1+idecal+nbold+i)
            end do
!
! CREATION DES NOUVEAUX OBJETS
            nbold = idecal+1+ndecal
            call wkvect(corres//'.PJEF_NB', 'V V I', nbnmes, iaconb)
            call wkvect(corres//'.PJEF_NU', 'V V I', nbold, iaconu)
            call wkvect(corres//'.PJEF_CF', 'V V R', nbold, iacocf)
! RECOPIE DES DONNEES
            do i = 1, nbnmes
                zi(iaconb-1+i) = zi(it1-1+i)
            end do
            do i = 1, nbold
                zi(iaconu-1+i) = zi(it2-1+i)
                zr(iacocf-1+i) = zr(it3-1+i)
            end do
!
! DESTRUCTION DES VECTEURS TAMPON
            call jedetr('TAMPON1')
            call jedetr('TAMPON2')
            call jedetr('TAMPON3')
!
!
            call jedetr('&&PJEFTE.LINONU1')
            call jedetr('&&PJEFTE.LINONU2')
!
        end do
!
    end if
!
!     AFFICHAGE DE LA CORRESPONDANCE DES NOEUDS
    ifres = iunifi('MESSAGE')
    write (ifres, '(A)') 'TABLE DE CORRESPONDANCE DES NOEUDS '
    write (ifres, '(A)') '----------------------------------'
    call jeveuo(corres//'.PJEF_NB', 'L', iaconb)
    call jeveuo(corres//'.PJEF_NU', 'L', iaconu)
    call jeveuo(corres//'.PJEF_CF', 'L', iacocf)
!
    kk = 0
    do i = 1, nbnmes
        labk8 = int_to_char8(i)
        write (ifres, 100) labk8
        nbno1 = zi(iaconb-1+i)
        do iocc = 1, nbno1
            kk = kk+1
            labk8 = int_to_char8(zi(iaconu-1+kk))
            coef = zr(iacocf-1+kk)
            write (ifres, 101) labk8, coef
        end do
    end do
!
100 format(' NOEUD MESURE :  ', a8)
101 format('       ', a8, '    POIDS : ', d12.5)
!
    call jedema()
end subroutine
