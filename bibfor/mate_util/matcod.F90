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
! person_in_charge: j-pierre.lefebvre at edf.fr
!
subroutine matcod(chmat, indmat, nbmat, imate, igrp, &
                  basename, codi, l_ther, base_)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/isnnem.h"
#include "asterf_types.h"
#include "asterfort/alfint.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupc.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveut.h"
#include "asterfort/jexatr.h"
#include "asterfort/tbexlr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: chmat, basename
    character(len=19) :: codi
    integer(kind=8) :: indmat, nbmat, imate, igrp
    aster_logical, intent(in) :: l_ther
    character(len=1), intent(in), optional :: base_
!
!-----------------------------------------------------------------------
!     MATERIAU CODE APPELE PAR RCMACO ET PMMACO
!-----------------------------------------------------------------------
!
!     BUT: CREER L'OBJET BASENAME//'      .CODI' ,LE REMPLIR ET RENVOYER
!          SON ADRESSE PAR RAPPORT A ZI
!
! IN  CHMAT  : NOM DU CHAM_MATER POUR ALFINT
! IN  INDMAT : INDICE DU PREMIER MATERIAU DANS LA LISTE
! IN  BASENAME : NOM SUR LEQUEL CONSTRUIRE LE CODI
! IN  NBMAT  : NOMBRE DE MATERIAUX DU NUMERO D'OCCURRENCE DE AFFE
! IN  IMATE  : NUMERO DU GROUPE
! IN  IGRP   : ADRESSE DU VECTEUR K8 CONTENANT LES NOMS DES MATERIAUX
! OUT  CODI   : OBJET MATERIAU CODE
!    CODI(1)   : ADRESSE ZK32  DE '.MATERIAU.NOMRC'
!    CODI(2)   : NOMBRE DE TYPES DE COMPORTEMENT
!
!         P.I = CODI(2+I)
!
!    CODI(2+I)  :POINTEUR DANS .CODI DU IEME COMPORTEMENT
!    CODI(P.I+0):NOMBRE DE COEFFICIENTS REELS
!    CODI(P.I+1):NOMBRE DE COEFFICIENTS COMPLEXES
!    CODI(P.I+2):NOMBRE DE COEFFICIENTS FONCTIONS
!    CODI(P.I+3):ADRESSE ZK16 RELATIVE AU .VALK (noms des parametres)
!    CODI(P.I+4):ADRESSE ZR  RELATIVE AU .VALR DES REELS
!    CODI(P.I+5):ADRESSE ZC  RELATIVE AU .VALC DES COMPLEXES
!    CODI(P.I+6):ADRESSE ZK16 RELATIVE AU .ORDR  (ou 1 si absent)
!    CODI(P.I+7):ADRESSE ZI   RELATIVE AU .KORD  (ou 1 si absent)
!
!         P.IF = P.I+LMAT
!
!    CODI(P.IF+LFCT*(K-1))  :NOMBRE DE POINTS DE LA FONCTION ASSOCIEE
!    CODI(P.IF+LFCT*(K-1)+1):ADRESSE ZK16 DU .PROL
!    CODI(P.IF+LFCT*(K-1)+2):ADRESSE ZR  DU .VALE
!    CODI(P.IF+LFCT*(K-1)+3):ADRESSE ZI  DU POINTEUR DE LONGUEUR(NAPPE)
!    CODI(P.IF+LFCT*(K-1)+4):ADRESSE ZR  DU .PARA (NAPPE)
!    CODI(P.IF+LFCT*(K-1)+5):LONUTI  DU .PARA (NAPPE)
!    CODI(P.IF+LFCT*(K-1)+6):POINTEUR SUPPLEMENTAIRE (TRACTION,TRC)
!    CODI(P.IF+LFCT*(K-1)+7):INDICE DE L'INTERVALLE POUR INTERPOLATION
!    CODI(P.IF+LFCT*(K-1)+8):INDICE SUPPLEMENTAIRE
!    CODI(P.IF+LFCT*(K-1)+9):coco = 1/2/3/4 : "code" du type de concept
!                            1 : fonction ou nappe
!                            2 : table TRC
!                            3 : liste de reels
!                            4 : liste de fonctions
!
!         P.IFC = CODI(P.IF+LFCT*(K-1)+6))
!
!    CODI(P.IFC)  :ADRESSE ZK8 DU .&&RDEP.PROL
!    CODI(P.IFC+1):ADRESSE ZR  DU .&&RDEP.VALE
!
! ----------------------------------------------------------------------
!
!
!
    integer(kind=8) :: iret, iretf, irett, iretlr, iretlf, nbcm, jnomrc, lmat, lfct, lsup
    integer(kind=8) :: jnbcm, l, nbv, nbtt, nbcot, nbcmt, nbco, iexi1, iexi2
    integer(kind=8) :: nbt, jdim, k, kk, nbk, lgcodi, isundf, idma, kr, nbreel
    integer(kind=8) :: imat, ipi, ipif, nbpts, ipifc, m, iretc
    integer(kind=8) :: jlisvr8, jlisvfo, jlisvi, jlisvr, code, iexi
    integer(kind=8) :: jcodi, jnomr, jjdim, jlcod, ipi0, nbfonc, kfonc, ipif2
    real(kind=8) :: tdef, prec
    character(len=4) :: knuma1
    character(len=3) :: knuma2
    character(len=3) :: knuma3
    character(len=6) :: k6
    character(len=8) :: nopara, nommat
    character(len=19) :: ch19, chma, listr, fon19
    character(len=1) :: base
! ----------------------------------------------------------------------
! PARAMETER ASSOCIE AU MATERIAU CODE
!
! --- LMAT   : NOMBRE DE PARAMETRES ASSOCIES AU COMPORTEMENT
! --- LFCT   : NOMBRE DE PARAMETRES ASSOCIES AUX FONCTIONS
! --- LSUP   : NOMBRE DE PARAMETRES SUPPLEMENTAIRE (COURBE &&RDEP)
    parameter(lmat=9, lfct=10, lsup=2)
! ----------------------------------------------------------------------
!
    call jemarq()

    if (present(base_)) then
        base = base_
    else
        base = 'V'
    end if
!
    call codent(imate, 'D0', knuma1)
    codi = ' '
    codi(1:8) = basename
    codi(9:13) = '.'//knuma1
    call jeexin(codi//'.CODI', iretc)
    if (iretc .ne. 0) then
        call jeveut(codi//'.CODI', 'L', jcodi)
        goto 999
    end if
!
    call wkvect('&&RCMACO.NBCM', 'V V I', nbmat, jnbcm)
    call wkvect('&&RCMACO.NOMR', 'V V I', nbmat, jnomr)
    call wkvect('&&RCMACO.JDIM', 'V V I', nbmat, jjdim)
    call wkvect('&&RCMACO.LCOD', 'V V I', nbmat, jlcod)
    do l = 1, nbmat
        nommat = zk8(igrp+indmat+l-1)
        call jeexin(nommat//'.MATERIAU.NOMRC', iret)
        ASSERT(iret .ne. 0)
        call jelira(nommat//'.MATERIAU.NOMRC', 'LONUTI', zi(jnbcm+l-1))
        call jeveut(nommat//'.MATERIAU.NOMRC', 'L', zi(jnomr+l-1))
        nbv = 0
        if (zk32(zi(jnomr+l-1)) .eq. 'ELAS_COQMU') nbv = 1
        if (zk32(zi(jnomr+l-1)+nbv) .eq. 'THER_COQMU') nbv = nbv+1
        if (nbv .gt. 0) zi(jnbcm+l-1) = nbv
    end do
!
    nbtt = 0
    nbcot = 0
    nbcmt = 0
    do l = 1, nbmat
        nbco = 0
        nbt = 0
        nommat = zk8(igrp+indmat+l-1)
        nbcm = zi(jnbcm+l-1)
        jnomrc = zi(jnomr+l-1)
        call codent(l, 'D0', knuma3)
        call jedetr('&&RCMACO.DIM'//knuma3)
        call wkvect('&&RCMACO.DIM'//knuma3, 'V V I', lmat*nbcm, zi(jjdim+l-1))
        jdim = zi(jjdim+l-1)
        do k = 1, nbcm
            kk = jdim+lmat*(k-1)
            call codent(k, 'D0', k6)
            chma = nommat//'.CPT.'//k6
            call codent(k, 'D0', knuma2)
            ch19 = chma(1:8)//'.'//knuma2//knuma1//knuma3
            call jedupc(' ', chma, 1, base, ch19, .false._1)
            call jelira(ch19//'.VALR', 'LONUTI', zi(kk))
            call jeveut(ch19//'.VALR', 'L', zi(kk+1))
            call jelira(ch19//'.VALC', 'LONUTI', zi(kk+2))
            call jeveut(ch19//'.VALC', 'L', zi(kk+3))
            call jelira(ch19//'.VALK', 'LONUTI', nbk)
            call jeveut(ch19//'.VALK', 'L', zi(kk+5))
            zi(kk+4) = (nbk-zi(kk)-zi(kk+2))/2
            nbco = nbco+zi(kk+4)
            if ((zk32(jnomrc+k-1) (1:8) .eq. 'TRACTION')) then
                zi(kk+6) = 1
                nbt = nbt+1
            end if
            if (zk32(jnomrc+k-1) (1:13) .eq. 'META_TRACTION') then
                zi(kk+6) = 1
                nbt = nbt+nbk/2
            end if
            call jeexin(ch19//'.ORDR', iexi)
            if (iexi .gt. 0) then
                call jeveut(ch19//'.ORDR', 'L', zi(kk+7))
                call jeveut(ch19//'.KORD', 'L', zi(kk+8))
            else
!               -- une adresse jeveux n'est jamais 1 :
                zi(kk+7) = 1
                zi(kk+8) = 1
            end if
        end do
        zi(jlcod+l-1) = 2+lmat*nbcm+lfct*nbco+lsup*nbt
        nbcmt = nbcmt+nbcm
        nbtt = nbtt+nbt
        nbcot = nbcot+nbco
    end do
!
    lgcodi = 2*nbmat+1+2*nbmat+lmat*nbcmt+lfct*nbcot+lsup*nbtt
    call wkvect(codi//'.CODI', base//' V I', lgcodi, jcodi)
    call jeveut(codi//'.CODI', 'E', jcodi)
    isundf = isnnem()
    do k = 1, lgcodi
        zi(jcodi+k-1) = isundf
    end do

    zi(jcodi) = nbmat

    ipi0 = 2*nbmat+1
    idma = ipi0
    do imat = 1, nbmat
        nommat = zk8(igrp+indmat+imat-1)
        nbcm = zi(jnbcm+imat-1)
        jnomrc = zi(jnomr+imat-1)
        jdim = zi(jjdim+imat-1)
        zi(jcodi+imat) = igrp+indmat+imat-1
        zi(jcodi+imat+nbmat) = idma
!
        zi(jcodi+idma) = jnomrc
        nbcm = zi(jnbcm+imat-1)
        zi(jcodi+idma+1) = nbcm
        ipi = jcodi+idma+2+nbcm
!
        do k = 1, nbcm
!
            chma = nommat//'.'//zk32(jnomrc+k-1) (1:10)
!
            kk = jdim+lmat*(k-1)
            zi(jcodi+idma+1+k) = ipi
            zi(ipi) = zi(kk)
            zi(ipi+1) = zi(kk+2)
            zi(ipi+2) = zi(kk+4)
            zi(ipi+3) = zi(kk+5)
            zi(ipi+4) = zi(kk+1)
            zi(ipi+5) = zi(kk+3)
            zi(ipi+6) = zi(kk+7)
            zi(ipi+7) = zi(kk+8)
!           zi(ipi+8) = zi(kk+8)
            ipif = ipi+lmat-1
!
! ---       boucle sur les coefficients reels :
!           ------------------------------------------
            do l = 0, zi(kk)-1
                ch19 = zk16(zi(kk+5)+l)
                if (ch19 .eq. 'PRECISION') prec = zr(zi(kk+1)+l)
            end do

            do l = 0, zi(kk)-1
                ch19 = zk16(zi(kk+5)+l)
                if (ch19 .eq. 'TEMP_DEF_ALPHA') then
                    tdef = zr(zi(kk+1)+l)
!
!                   boucle sur les fonctions :
!                   ------------------------
                    do m = 0, zi(kk+4)-1
                        ch19 = zk16(zi(kk+5)+zi(kk)+zi(kk+2)+zi(kk+4)+m)
                        nopara = zk16(zi(kk+5)+zi(kk)+zi(kk+2)+m) (1:8)
                        if (nopara(1:5) .eq. 'ALPHA' .or. nopara .eq. 'F_ALPHA ' .or. &
                            nopara .eq. 'C_ALPHA ') then
!
!                           interpolation des coefficients de dilatation alpha
!                           en tenant compte de la temperature de definition tdef :
!                           -----------------------------------------------------
                            if (chmat .ne. '&chpoint') then
                                call alfint(chmat, imate, nommat, tdef, nopara, &
                                            k, prec, ch19, l_ther)
                                zk16(zi(kk+5)+zi(kk)+zi(kk+2)+zi(kk+4)+m) = ch19(1:16)
                            end if
                        end if
                    end do
                end if
            end do

!          -- boucle sur les coefficients fonction/table/lisv :
!          -----------------------------------------------------------
            do l = 0, zi(kk+4)-1
                ch19 = zk16(zi(kk+5)+zi(kk)+zi(kk+2)+zi(kk+4)+l)
                call exisd('FONCTION', ch19(1:8), iretf)
                call exisd('TABLE', ch19(1:8), irett)

!               -- cas des LISV_R8 / LISV_FO :
                iretlr = 0
                iretlf = 0
                call jeexin(ch19(1:16)//'.LISV_R8', iexi1)
                call jeexin(ch19(1:16)//'.LISV_FO', iexi2)
                if (iexi1 .gt. 0) iretlr = 1
                if (iexi2 .gt. 0) iretlf = 1

                ASSERT(iretf .eq. 1 .or. irett .eq. 1 .or. iretlr .eq. 1 .or. iretlf .eq. 1)

!               -- des fonctions sont creees sur la volatile (routine alfint) ---
                if (iretf .eq. 1) then
!                   -- cas des fonctions :
                    call jeveut(ch19//'.PROL', 'L', zi(ipif+1))
                    zi(ipif+7) = 1
                    zi(ipif+8) = 1
                    zi(ipif+9) = 1
                    if (zk24(zi(ipif+1)) (1:1) .eq. 'C' .or. zk24(zi(ipif+1)) (1:1) .eq. 'F') then
                        call jeveut(ch19//'.VALE', 'L', zi(ipif+2))
                        call jelira(ch19//'.VALE', 'LONMAX', nbpts)
                        zi(ipif) = nbpts/2
                    else if (zk24(zi(ipif+1)) (1:1) .eq. 'N') then
                        call jeveut(ch19//'.VALE', 'L', zi(ipif+2))
                        call jeveut(jexatr(ch19//'.VALE', 'LONCUM'), 'L', zi(ipif+3))
                        call jeveut(ch19//'.PARA', 'L', zi(ipif+4))
                        call jelira(ch19//'.PARA', 'LONUTI', zi(ipif+5))
                    else if (zk24(zi(ipif+1)) (1:1) .eq. 'I') then
                    else
                        call utmess('F', 'MODELISA6_64', sk=zk24(zi(ipif+1)))
                    end if

                else if (irett .eq. 1) then
!                   -- cas des tables (TRC) :
                    listr = '&&'//ch19(1:8)//'_LR8'
                    call jeexin(listr//'.VALE', iretc)
                    if (iretc .eq. 0) then
                        call tbexlr(ch19, listr, 'V')
                    end if
                    call jeveut(listr//'.VALE', 'L', zi(ipif))
                    zi(ipif+1) = 0
                    zi(ipif+2) = 0
                    zi(ipif+9) = 2

                else if (iretlr .eq. 1) then
!                   -- cas des mater_LISV_R8 :
                    code = -1
                    call jeveuo(ch19(1:16)//'.LISV_R8', 'L', jlisvr8)
                    call jelira(ch19(1:16)//'.LISV_R8', 'LONMAX', nbreel)
                    call jedetr(ch19(1:16)//'.LISV_VR')
                    call wkvect(ch19(1:16)//'.LISV_VR', 'V V R', nbreel+1, jlisvr)
                    call jeveut(ch19(1:16)//'.LISV_VR', 'E', jlisvr)
                    zr(jlisvr-1+1) = dble(nbreel)
                    call jedetr(ch19(1:16)//'.LISV_IA')
                    call wkvect(ch19(1:16)//'.LISV_IA', 'V V I', 2, jlisvi)
                    call jeveut(ch19(1:16)//'.LISV_IA', 'E', jlisvi)
                    zi(jlisvi-1+1) = code
                    zi(jlisvi-1+2) = jlisvr
                    do kr = 1, nbreel
                        zr(jlisvr-1+1+kr) = zr(jlisvr8-1+kr)
                    end do

                    zi(ipif) = jlisvi
                    zi(ipif+1) = 0
                    zi(ipif+2) = 0
                    zi(ipif+9) = 3

                else if (iretlf .eq. 1) then
!                   -- cas des mater_LISV_FO :
                    code = -2
!                   -- on cree un vecteur d'entiers pour stocker les adresses et les infos
!                      necessaires a l'evaluation rapide des listes de fonctions :
                    call jeveuo(ch19(1:16)//'.LISV_FO', 'L', jlisvfo)
                    call jelira(ch19(1:16)//'.LISV_FO', 'LONMAX', nbfonc)
                    call jedetr(ch19(1:16)//'.LISV_VR')
                    call wkvect(ch19(1:16)//'.LISV_VR', 'V V R', nbfonc+1, jlisvr)
                    call jeveut(ch19(1:16)//'.LISV_VR', 'E', jlisvr)
                    zr(jlisvr-1+1) = dble(nbfonc)
                    call jedetr(ch19(1:16)//'.LISV_IA')
                    call wkvect(ch19(1:16)//'.LISV_IA', 'V V I', 3+lfct*nbfonc, jlisvi)
                    call jeveut(ch19(1:16)//'.LISV_IA', 'E', jlisvi)
                    zi(jlisvi-1+1) = code
                    zi(jlisvi-1+2) = jlisvr
                    zi(jlisvi-1+3) = nbfonc
                    do kfonc = 1, nbfonc
                        fon19 = zk8(jlisvfo-1+kfonc)
                        ipif2 = jlisvi+3+lfct*(kfonc-1)
                        call jeveut(fon19//'.PROL', 'L', zi(ipif2+1))
                        zi(ipif2+7) = 1
                        zi(ipif2+8) = 1
                        if (zk24(zi(ipif2+1)) (1:1) .eq. 'C' &
                            .or. zk24(zi(ipif2+1)) (1:1) .eq. 'F') then
                            call jeveut(fon19//'.VALE', 'L', zi(ipif2+2))
                            call jelira(fon19//'.VALE', 'LONMAX', nbpts)
                            zi(ipif2) = nbpts/2
                        else if (zk24(zi(ipif2+1)) (1:1) .eq. 'N') then
                            call jeveut(fon19//'.VALE', 'L', zi(ipif2+2))
                            call jeveut(jexatr(fon19//'.VALE', 'LONCUM'), 'L', zi(ipif2+3))
                            call jeveut(fon19//'.PARA', 'L', zi(ipif2+4))
                            call jelira(fon19//'.PARA', 'LONUTI', zi(ipif2+5))
                        else if (zk24(zi(ipif2+1)) (1:1) .eq. 'I') then
                        else
                            call utmess('F', 'MODELISA6_64', sk=zk24(zi(ipif2+1)))
                        end if
                    end do

                    zi(ipif) = jlisvi
                    zi(ipif+1) = 0
                    zi(ipif+2) = 0
                    zi(ipif+9) = 4

                else
                    call utmess('F', 'MODELISA6_64', sk=ch19(1:8))
                end if
!
                if (zi(kk+6) .eq. 1) then
                    if ((zk32(jnomrc+k-1) (1:8) .eq. 'TRACTION') .or. &
                        (zk32(jnomrc+k-1) (1:13) .eq. 'META_TRACTION')) then
                        ipifc = ipif+lfct
                        zi(ipif+6) = ipifc
                        ch19 = nommat//'.&&RDEP'
                        call jeveut(ch19//'.PROL', 'E', zi(ipifc))
                        call jeveut(ch19//'.VALE', 'E', zi(ipifc+1))
                        ipif = ipifc+lsup
                    else
                        ipif = ipif+lfct
                    end if
                else
                    ipif = ipif+lfct
                end if
            end do
            ipi = ipif
        end do
        idma = idma+zi(jlcod+imat-1)
    end do
!
! --- MENAGE
    do l = 1, nbmat
        call codent(l, 'D0', knuma3)
        call jedetr('&&RCMACO.DIM'//knuma3)
    end do
    call jedetr('&&RCMACO.NBCM')
    call jedetr('&&RCMACO.NOMR')
    call jedetr('&&RCMACO.JDIM')
    call jedetr('&&RCMACO.LCOD')
!
!
999 continue
    call jedema()
!
end subroutine
