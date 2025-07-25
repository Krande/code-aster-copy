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

subroutine nmexso(mesh, ds_inout, sddyna, nume_dof)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/ndynkk.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_InOut), intent(in) :: ds_inout
    character(len=19), intent(in) :: sddyna
    character(len=24), intent(in) :: nume_dof
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Computation
!
! Initialization of FORCE_SOL load
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ds_inout         : datastructure for input/output management
! In  sddyna           : name of dynamic parameters datastructure
! In  nume_dof         : name of numbering object (NUME_DDL)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=15) :: sdexso
    character(len=19) :: sdexsz
    character(len=8) :: cnfsol, result, nomacr
    character(len=24) :: magrno, maille, nprno
    character(len=24) :: tabequ, tabfrq, tabinf, nomres
    integer(kind=8) :: ieqint, jfrq, iddint, jnomre
    character(len=24) :: gnintf, tabrig, tabmas, tabamo
    integer(kind=8) :: jrig, jmas, jamo
    integer(kind=8) :: gd, aprno
    real(kind=8) :: pasa, pasm, pas, ainst, rinst, nbm
    integer(kind=8) :: idno
    integer(kind=8) :: ibid
    integer(kind=8) :: ifreq, i1, i2, inoe, ino, ima, iddl, icmp
    character(len=24) :: uniamo, unirig, unimas, unifor, ipasm
    integer(kind=8) :: unitea, uniter, unitem, unitef, npasm, ibin
    integer(kind=8) :: unifrq
    integer(kind=8) :: nfreq, nfreqm, nfreqa
    integer(kind=8) :: nbmode, nbmod2, nbno, nddint, ncmp, nec, neq
    character(len=24) :: nchsol
    integer(kind=8) :: jchsol
    character(len=8), pointer :: vnomacr(:) => null()
    character(len=24), pointer :: veiss(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Get parameters
!
    result = ds_inout%result
!
! --- INITIALISATIONS
!
    tabfrq = '&&NMEXSO.FREQ'
!
! --- ACCES NUMEROTATION
!
    call dismoi('NB_EQUA', nume_dof, 'NUME_DDL', repi=neq)
    call dismoi('NUM_GD_SI', nume_dof, 'NUME_DDL', repi=gd)
    nec = nbec(gd)
!
! --- ACCES SD EXCIT_SOL
!
    call ndynkk(sddyna, 'SDEXSO', sdexsz)
    sdexso = sdexsz(1:15)
!
! --- SAUVEGARDE NOM DU RESULTAT
!
    nomres = sdexso(1:15)//'.RESU'
    call wkvect(nomres, 'V V K8', 1, jnomre)
    zk8(jnomre) = result
!
! --- RECUPERATION CHARGE
!
    nchsol = sdexso(1:15)//'.CHAR'
    call jeveuo(nchsol, 'L', jchsol)
    cnfsol = zk8(jchsol)
    call jeveuo(cnfsol//'.CHME.VEISS', 'L', vk24=veiss)
!
! --- NOMS FICHIERS
!
    unirig = veiss(1)
    unimas = veiss(2)
    uniamo = veiss(3)
    unifor = veiss(4)
    uniter = 0
    unitem = 0
    unitea = 0
    unitef = 0
!
! --- GROUP_NO_INTERF
!
    gnintf = veiss(5)
    if (gnintf .eq. ' ') then
        maille = veiss(6)
        call jeveuo(mesh//'.NOMACR', 'L', vk8=vnomacr)
        call jenonu(jexnom(mesh//'.SUPMAIL', maille), ima)
        nomacr = vnomacr(ima)
        call jelira(nomacr//'.LINO', 'LONMAX', ival=nbno)
        call jeveuo(nomacr//'.LINO', 'L', idno)
    else
        magrno = mesh(1:8)//'.GROUPENO'
        call jelira(jexnom(magrno, gnintf), 'LONUTI', ival=nbno)
        call jeveuo(jexnom(magrno, gnintf), 'L', idno)
    end if
!
! --- NOMBRE DE DDL INTERNES
!
    nprno = nume_dof(1:14)//'.NUME.PRNO'
    call jenonu(jexnom(nprno(1:19)//'.LILI', '&MAILLA'), ibid)
    call jeveuo(jexnum(nprno, ibid), 'L', aprno)
    nddint = 0
    do ino = 1, nbno
        inoe = zi(idno+ino-1)
        ncmp = zi(aprno+(nec+2)*(inoe-1)+2-1)
        nddint = nddint+ncmp
    end do
!
! --- TABLEAU DES NUMEROS D EQUATION ACTIFS DE L INTERFACE
!
    tabequ = sdexso(1:15)//'.EQINT'
    call wkvect(tabequ, 'V V I', nddint, ieqint)
    iddint = 0
    do ino = 1, nbno
        inoe = zi(idno+ino-1)
        ncmp = zi(aprno+(nec+2)*(inoe-1)+2-1)
        iddl = zi(aprno+(nec+2)*(inoe-1)+1-1)
        do icmp = 1, ncmp
            iddint = iddint+1
            zi(ieqint+iddint-1) = iddl+icmp-1
        end do
    end do
!
! --- OUVERTURE DES FICHIERS
!
    if (unirig .ne. ' ') read (unirig, '(I24)') uniter
    if (unimas .ne. ' ') read (unimas, '(I24)') unitem
    if (uniamo .ne. ' ') read (uniamo, '(I24)') unitea
    if (unifor .ne. ' ') read (unifor, '(I24)') unitef
!
! --- QUEL FICHIER VA DONNER LA FREQUENCE ?
!
    unifrq = 0
    if (uniter .eq. 0) then
        if (unitem .eq. 0) then
            unifrq = unitea
        else
            unifrq = unitem
        end if
    else
        unifrq = uniter
    end if
!
    if (unifrq .eq. 0) then
        ASSERT(.false.)
    end if
    if (unifrq .eq. uniter) then
        call utmess('I', 'DYNAMIQUE_20')
    end if
    if (unifrq .eq. unitem) then
        call utmess('I', 'DYNAMIQUE_21')
    end if
    if (unifrq .eq. unitea) then
        call utmess('I', 'DYNAMIQUE_22')
    end if
!
! --- FICHIER BINAIRE ?
    read (veiss(8), '(I24)') ibin
!
! --- LECTURE DU PAS D'ACTUALISATION
!
    if (ibin .eq. 0) then
        rewind unifrq
        read (unifrq, *) ainst, pas
    else
        open (unit=unifrq, form='unformatted', status='old', access='stream')
        read (unifrq) ainst, pas, nbm
    end if
    nfreq = nint(ainst)
    call utmess('I', 'DYNAMIQUE_23', si=nfreq, sr=pas)
!
! --- TABLEAU DES FREQUENCES
!
    call wkvect(tabfrq, 'V V R', nfreq, jfrq)
    if (ibin .eq. 0) then
        do ifreq = 1, nfreq
            zr(jfrq+ifreq-1) = (ifreq-1)*pas
        end do
    else
        read (unifrq) (zr(jfrq+ifreq-1), ifreq=1, nfreq)
    end if
!
! --- VERIFICATIONS
!
    if (unitem .ne. 0) then
        if (unifrq .ne. unitem) then
            if (ibin .eq. 0) then
                rewind unitem
                read (unitem, *) ainst, pasm
            else
                open (unit=unitem, form='unformatted', status='old', access='stream')
                read (unitem) ainst, pasm, nbm
            end if
            nfreqm = int(ainst)
            if (nfreqm .ne. nfreq) then
                call utmess('F', 'DYNAMIQUE_30')
            end if
            if (ibin .ne. 0) read (unitem) (zr(jfrq+ifreq-1), ifreq=1, nfreq)
        end if
    end if
    if (unitea .ne. 0) then
        if (unifrq .ne. unitea) then
            if (ibin .eq. 0) then
                rewind unitea
                read (unitea, *) ainst, pasa
            else
                open (unit=unitea, form='unformatted', status='old', access='stream')
                read (unitea) ainst, pasa, nbm
            end if
            nfreqa = int(ainst)
            if (nfreqa .ne. nfreq) then
                call utmess('F', 'DYNAMIQUE_31')
            end if
            if (ibin .ne. 0) read (unitea) (zr(jfrq+ifreq-1), ifreq=1, nfreq)
        end if
    end if
!
! --- TAILLE TABLEAUX
!
    nbmode = nddint
    nbmod2 = nbmode*nbmode
!
! --- SAUVEGARDE INFORMATIONS
!
    tabinf = sdexso(1:15)//'.TABI'
    call wkvect(tabinf, 'V V R', 4, iddint)
    zr(iddint-1+1) = pas
    zr(iddint-1+2) = unitef
    zr(iddint-1+3) = nddint
    ipasm = veiss(7)
    read (ipasm, '(I24)') npasm
    zr(iddint-1+4) = npasm
!
! --- LECTURE MATRICE REDUITE RIGIDITE A L'INTERFACE
!
    tabrig = sdexso(1:15)//'.RIGT'
    call wkvect(tabrig, 'V V R', nbmod2*nfreq, jrig)
    if (uniter .ne. 0) then
        if (ibin .eq. 0) then
            do ifreq = 1, nfreq
                read (uniter, *) rinst
                read (uniter, 100) &
                    ((zr(jrig+(ifreq-1)*nbmod2+(i2-1)*nbmode+i1-1), i2=1, nbmode), i1=1, nbmode)
            end do
        else
            do ifreq = 1, nfreq
                read (uniter) &
                    ((zr(jrig+(ifreq-1)*nbmod2+(i2-1)*nbmode+i1-1), i2=1, nbmode), i1=1, nbmode)
            end do
            close (unit=uniter)
        end if
    end if
!
! --- LECTURE MATRICE REDUITE MASSE A L'INTERFACE
!
    tabmas = sdexso(1:15)//'.MAST'
    call wkvect(tabmas, 'V V R', nbmod2*nfreq, jmas)
    if (unitem .ne. 0) then
        if (ibin .eq. 0) then
            do ifreq = 1, nfreq
                read (unitem, *) rinst
                read (unitem, 100) &
                    ((zr(jmas+(ifreq-1)*nbmod2+(i2-1)*nbmode+i1-1), i2=1, nbmode), i1=1, nbmode)
            end do
        else
            do ifreq = 1, nfreq
                read (unitem) &
                    ((zr(jmas+(ifreq-1)*nbmod2+(i2-1)*nbmode+i1-1), i2=1, nbmode), i1=1, nbmode)
            end do
            close (unit=unitem)
        end if
    end if
!
! --- LECTURE MATRICE REDUITE AMORTISSEMENT A L'INTERFACE
!
    tabamo = sdexso(1:15)//'.AMOT'
    call wkvect(tabamo, 'V V R', nbmod2*nfreq, jamo)
    if (unitea .ne. 0) then
        if (ibin .eq. 0) then
            do ifreq = 1, nfreq
                read (unitea, *) rinst
                read (unitea, 100) &
                    ((zr(jamo+(ifreq-1)*nbmod2+(i2-1)*nbmode+i1-1), i2=1, nbmode), i1=1, nbmode)
            end do
        else
            do ifreq = 1, nfreq
                read (unitea) &
                    ((zr(jamo+(ifreq-1)*nbmod2+(i2-1)*nbmode+i1-1), i2=1, nbmode), i1=1, nbmode)
            end do
            close (unit=unitea)
        end if
    end if
!
    call jedetr(tabfrq)
100 format((6(1x, 1pe13.6)))
    call jedema()
end subroutine
