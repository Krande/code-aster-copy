# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

import os.path as osp
import unittest

from code_aster.MacroCommands.Miss.miss_fichier_cmde import MissCmdeGen, remove_comments
from code_aster.MacroCommands.Miss.miss_resu_aster import ResuAsterReader
from code_aster.MacroCommands.Miss.miss_resu_miss import MissCsolReader
from code_aster.MacroCommands.Utils import test_utils
from code_aster.Objects.table_py import Table


class TestMissCmde(unittest.TestCase):

    """test generator of miss commands file"""

    def setUp(self):
        """set up parameters"""

        class Parameters(dict):

            """fake MISS_PARAMETERS for unittests"""

            def set(self, key, value):
                """assign a value"""
                self[key] = value

        class Struct:

            """fake structure"""

            titre = "PRODUIT PAR CALC_MISS"

        self.struct = Struct()
        self.par = Parameters()
        self.par.update(
            {
                # valeurs par défaut des mots-clés, cf. calc_miss.capy
                "PROJET": "MODELE",
                "FREQ_MIN": None,
                "FREQ_MAX": None,
                "FREQ_PAS": None,
                "LIST_FREQ": None,
                "FREQ_IMAG": None,
                "Z0": 0.0,
                "SURF": "NON",
                "ISSF": "NON",
                "ALLU": 0.0,
                "RFIC": 0.0,
                "ALGO": None,
                "DREF": None,
                "SPEC_MAX": None,
                "SPEC_NB": None,
                "OFFSET_MAX": None,
                "OFFSET_NB": None,
                "TYPE": "ASCII",
                "AUTO": "NON",
                "OPTION_DREF": "NON",
                "OPTION_RFIC": "NON",
                "COEF_OFFSET": 12,
                "_hasPC": False,
                "_hasSL": False,
                "_nbPC": 0,
                "_auto_first_LT": None,
            }
        )
        self._debug = False

    def fname(self, ext):
        """use PROJET"""
        return "./" + self.par["PROJET"] + "." + ext

    def test01_sdlx103a(self):
        """use LIST_FREQ"""
        refe = """
eJzdVsuOmzAU3fsr7jpVUsho2m4pOAwqr4IZjboZUeIklnikPKr+Ur+jP9ZrJ9OQCQm0yyKi2LJ9
7rnnPswMLP61a2HNIQfedvjviaYhM5gPPWSGK35VwJaXvBbfOo4nG9iIbCd43YDnxPG1s30Ym/o0
gthyn3TtzlCoRx4dFLzswDKYMQhEDiu4xkRb8x5xMrvczBwWURJGgZU4DEIjAtNwzecDTcRwedZ2
EgXNpiLP0+0QjETyHgMXFm9fKC+K71V+5L0RpWhFVSoptnXV7flVAdH1KEhCAgBLQMjEo2Tl+PJn
D8EV1foqGDEfDC8krhPRPrNst1dIYVqnBUeNGliDKFu+rVOJeyM+xPEZtSODOYEPETUZvIMPwCLH
8G2Xgr7E95Il5CnspXRyvKk5ZkWZ3ZDglc1VRD8n1DcpqOdOaqMvlst7TdPoG00/TOWsN33/sjok
WlN1zXxdYUTLURoX57uz49C0dadyZKo7sRV4huMrb3RQ0T6OMeTkk+nJWMejZmV2jRg8M7U8mZrr
JGaR0YqN4CdjK1GOVZfcbPWLUdXWD551N/JG7Td3ab3liK0q+G+cwILuOTGcXBiDVHqTpSetJiCf
DMRLD75o0sL9QtOOmYTzi9pphos626XFvsEiysQafZyQ2lhJpoO5bIWugb7Fj5i0AxSItgB89cWf
jQ8TN4Zj+1RY0jzr8n9wgdAnaoL04b8Mb0+XTVVm0mwjydg157f6Y0+bOFxFl1g19j28jUAUe75O
p3XBA15iMxscL6TqGKyCyOSnuh0vxPM7dF817RxlFa2K2q2Wr66K8+3KGXTh18+JPpAwiBl2D/Oh
rzbePN2zVALgo+MbBL2zTNXugQUJozFeYGEix2rEDiPsQeFVUkeFG9EU8uNj8h1znZ9EJFJtesbM
stwDncR5RevYR6cq/CqGI19Y6lvgN0yQXTY=
"""
        self.par.update(
            {
                "PROJET": "SDLX103A",
                "LIST_FREQ": (12.25, 12.50, 12.75),
                "TYPE": "BINAIRE",
                "Z0": 5.0,
                "DREF": 1.0,
                "ALGO": "REGU",
                "OFFSET_MAX": 20,
                "OFFSET_NB": 200,
            }
        )
        gen = MissCmdeGen(self.par, self.struct, self.fname)
        txt = gen.build()
        if self._debug:
            with open("/tmp/test01_sdlx103a.in", "w") as f:
                f.write(txt)
        diff = self._diffcompress(refe, txt)
        assert diff.strip() == "", diff

    def test02_zzzz200b(self):
        """use FREQ_MIN/FREQ_MAX/FREQ_PAS"""
        refe = """
eJzdVtmumzAQffdXzHOqpEmqVn1F4HBR2WrMVdWXK0qcxBJLylL1l/od/bGOTXTJQgLtYxERVhif
OXNmMTOwxLe2ga2ADETT4tOTdU1mMB+6yAzf+GUOe1GISn5vBe6sYSfTgxRVDZ4TRff2nsPY1KcM
vMCiLtWYJxYt5KJowTK4MQhDujf4jsumEme0yezWmDucURKywIodDqHBwDRc86UjiRiuSJtWoaDb
RGZZsh+CUUjec+DC4m1HeJH/KLMT650sZCPLQsuwr8r2KO6Kh2GzIA4JAKwBAWOPko3jq589BJeX
27tgxHwyvJC4DqM9r/Rw1DhhUiW5QH1q2IIsGrGvEoX6IDPE8Tm1mcGdwAdGTQ4f4CNw5hi+7VJY
rfG+5QhZAkclm1rvKoH1UKQPBLjyuWH0c0x9k4JFURRYLZb6om+WSzCu/llhAqNrqyHZ6rKt59sS
M1qMUrnZ315sh7qpWl0jU0OKrMAzHL8LB3S+T2tMOvlkeirb0ahbVV8jDi9crXtX8xWJODMauZOi
d7aRxVh3KWPrvBl1b/0UafugdrS9eUiqvUBs3cF/EwQ29FkQwwWGOUhUNGnSazUBuXcQrT34ulTS
fFgs378WGEQx21x1UD3c2OkhyY81tlIqtxjlhALHfjIdeEes0DUwuugZy3aABFkuAO/V4tXwaaJh
OGanE5NkaZv9QwiEfqEmqBj+0wSfKbMri1Q5rhUduxLi0Zw8UycKN+wWq8L5hycSyPwotsm0adjh
xTa3wfFCqrfBJmCm6Ht3vBkvz9FjWTdzFFY2Om+PRr8+Mi7NdTAYwu9fE2MgYRBxnCDmU681nj/t
i9IBgGBclqkHPvAg5jTCAyyM1VqveLfCCRTepXPStpZ1rj48Jp8y95gpPKJUphe8LMvtyMTOFanT
DJ2q7FXuRr6t9JfAH4AuXdY=
"""
        self.par.update(
            {
                "FREQ_MIN": 1.0,
                "FREQ_MAX": 10.0,
                "FREQ_PAS": 1.0,
                "Z0": -6.05,
                "SURF": "OUI",
                "ALGO": "DEPL",
                "OFFSET_MAX": 60,
                "OFFSET_NB": 400,
                "SPEC_MAX": 0.03,
                "SPEC_NB": 2048,
            }
        )
        gen = MissCmdeGen(self.par, self.struct, self.fname)
        txt = gen.build()
        if self._debug:
            with open("/tmp/test02_zzzz200b.in", "w") as f:
                f.write(txt)
        diff = self._diffcompress(refe, txt)
        assert diff.strip() == "", diff

    def test03_fdlv112b(self):
        """use ISSF"""
        refe = """
eJzdVt1umzAUvvdTnOtUyUi6TdslAZMiEWDGRNNuKkqcxBKBlJ9pr7Tn2Ivt2KQNaX5brTdDRXGN
fc53vvPbA1s8NDXMBWQg6gZ/p7KqSA/6xx7Swy9+sYalyEUpHxuBNytYyHQlRVnB1I2iU3e7YibU
pwwc25sNh6OxlrrF0cBa5A3YJjePCiLtF/zGZV2KDnDSOzzMXc4oCVlgxy6H0GRgmZ5138JEGZ5I
60ZJQbWJzLJkeUyMkjSdBR4MPjxBHqx/FtkW90LmspZFrqlYlkWzEScJRNNZEIfEcf3JsevrYn7y
MrHuzGlIPJfRLpJ0tdGSwqRM1gI5qWAOMq/FskyU3DP+IK7P6YSZ3A18YNTi8Bm+AGeu6U88CsMR
/h2ihCyBjaJKrRelwCjI0zMmv9DpMPotpr5FwaaAz3Bg6If2jSGYauf2aefGMNBpkdobDb4+7x0j
riqaqj8v0Iv5RSgH95u96/hPdtEYEtnB1HT91gTQXlXr/lAZgO9HEnFm1nIhhfJ29M+VjjpKR/je
EtvxYhdJZcEerzfGLVjUw4TTu592u2PKNeFGl/AoZg7u7EA7Mr+Ul+qw3U1jnZW/RNqciUB93lol
5VKgbJ37ryEDS0HHA8fDtKpVDmB5SnacXyF5pyAaTeGHoTR86nIEB1lYHS8H6SpZbypMx1TO0cYr
kgRz0nIBnRl6JtoWzdBlRyAQY4BOwm/PB++uPBheOqfdkmRpk73BBEK/UwuUDWfcu8gaqQr2azzc
FqJx0aQZxkdTQqYaz1XVB8U4TDwCD2JOI4jMGYXpzLFxg0PsAsfXnYY2OAGz6H8Zlh1/Loo8VWor
BWZSCnGuQ7QKtU+j0GGHskrkHvsvyPVGzJPr+kArL57wiaKd6muae/EuMfPeeJlAlpsnT27K4iFD
9NsmXCbphSZ8qGniBWPT29VeHTMPb4n8rmTfdvDOntQLJXp/LtsUVd3HwJW19s25sUKPI/vHNf1I
+p/fV7JOwiDi2Fesu24843TT3CvfYayrlLX0OLHNbByQwlit9Yq3K+xM4UlA23ioZLVWw+zVM8xp
bEoiaetIF5ltey0cVW72YG2767XsvvDfhYldqSF/Ac9u8B8=
"""
        self.par.update(
            {
                "PROJET": "FDLV112B",
                "FREQ_MIN": 0.1,
                "FREQ_MAX": 3.0,
                "FREQ_PAS": 2.9,
                "Z0": 5.0,
                "SURF": "NON",
                "ALGO": "REGU",
                "OFFSET_MAX": 1000,
                "OFFSET_NB": 5000,
                "ISSF": "OUI",
                "MATER_FLUIDE": {
                    "RHO": 1000.0,
                    "CELE": 1500,
                    "AMOR_BETA": 0.0,
                    "DEMI_ESPACE": "OUI",
                },
                "_hasSL": True,
            }
        )
        gen = MissCmdeGen(self.par, self.struct, self.fname)
        txt = gen.build()
        if self._debug:
            with open("/tmp/test03_fdlv112b.in", "w") as f:
                f.write(txt)
        diff = self._diffcompress(refe, txt)
        assert diff.strip() == "", diff

    def test04_miss03c_inci(self):
        """first file built when FICHIERS_TEMPS is enabled"""
        refe = """
eJzdVtuOmzAQffdXzHOq0CTVVvuKwGFRudWYVdWXFSVOYolLiqHq53fs3W1IwoZ0H4uIIDA+M+fM
xczAFT/6DjYCShBdj9dQKkVmMB87yAzfRE0FO1GLVv7sBa5UsJXFXopWQein6VtrhzAejSgznp6C
/FDmhTDIL7H0UIm6B9fm9igYeX6D77jsWjEInswujbnPGSUJi93M55DYDBw7cJ6eQ0WMQBRdr1HQ
bS7LMt+NwWik8DEOwPo4DNuqfjXlS+xbWctONrWRZNc2/UG8KSRKwOIsIQCwAoTNQkrWfqR/3hhc
1WzeBCPOgx0mJPAZPY+u2B8MWpK3eSVQKwUbkHUndm2usa/kivgRpx6zuR9HwKjD4TPcA2e+HXkB
heUKz8tIoczhoCXU99tWYIXUxRUZznyuGf2a0cih4FKUBu6shTnofLEEWz9ZvT75gE8SOz23GhNP
Nb2abxrMbj0ZysX6/mQ5qK7tTb3cSil149D2I0NnCSbrL/dIhnxxQp3zdNKtrrIJhyeuVkdX8yVJ
ObM7uZXi6Gwt66lO08busDFNn/0WRX+ldoy9s8/bnUBs083/QgKbe0BivMAwB7lmU+RHrW5APjpI
VyF8X2hp7qz713LC/2nG1qN9pMabvNjn1UFhQxVyg1xvKHPsKseHT8RNAhs5po+wtEZCIQsL8Fxa
fw0fbjRMpuxMevKy6Mt3UCD0G3VAc/iv0zzQZ9vUhXavdFBeK8S1mTnQKE3W7BKrLfQ4kqrSW+dN
Wmce92AdM0ccu3a6DU9300OjujmKKTuTq2tD32wWp+bvCP2Il8QpxzHiPJxLjVtR/6RRiSZHzeQH
HmecpsR1A32bksx/vuIcSoZD61ZCZ5JNfN6YDfgPEoFDsQ==
"""
        self.par.update(
            {
                "PROJET": "Miss_Laplace",
                "ALGO": "DEPL",
                "DREF": 5,
                "FREQ_IMAG": 1.0,
                "INST_FIN": 2.0,
                "OFFSET_MAX": 40,
                "OFFSET_NB": 400,
                "PAS_INST": 0.05,
                "RFIC": 0.0,
                "SPEC_MAX": 0.075,
                "SPEC_NB": 16384,
                "SURF": "OUI",
                "Z0": -5.8,
            }
        )
        gen = MissCmdeGen(self.par, self.struct, self.fname, lapl_temps=True)
        txt = gen.build()
        if self._debug:
            with open("/tmp/test04_miss03c_inci.in", "w") as f:
                f.write(txt)
        diff = self._diffcompress(refe, txt)
        assert diff.strip() == "", diff

    def test05_miss03c_loop(self):
        """using LAPL_TEMPS"""
        refe = """
eJzdVtuO2jAQffdXzDMVlMB2u32MEpONmlsTZ1X1ZZUGA5ZyoblU/aV+R3+sY4NKgEDSPjYCYeHx
mTNnLs4ETP61bWDNIQPetPjriromE5j2PWSCO16Zw5YXvBLfWo4na9iIdCd4VYNrR9Gts10Yi3o0
VJ5enWSfJSlXyEcuLeS8aMHUmd4LRg47uMdEU/EOeTK5NmY2CykJQt+MbQaBHoKhO8brgSpiODxt
WomCbhORZcm2D0YiuS++A7O3Xdqz/HuZHblvRCEaURZKkm1Vtnt+U0iUIPTjgADAAhA2dilZ2Z78
Wn1webm+CUaMZ90NiGOH9JJdutsrtCCpkpyjVjWsQRQN31aJxL6TK2J7jFqhzmzfg5AaDB7hCVho
657lUNAW+LlmClkCeymhXG8qjhVSpHdkuPC5CumnmHoGBfVoUh9tNl8sl09L+mY+J7Zr+XLncfbw
4f3D4zv8U+vTqy7berouMaHFoPer8+3ZcaibqlUlMjaKyPRd3fZUEBqoRB/XmG3y0XBlmqNBt7Kw
BhyeuVqcXE01ErFQb8RG8JOzlSiGmksam91eVK31g6ftnXJR9sYuqbYcsVUD/00Q2M+dIPprCnOQ
yGjS5KTVCOSTg2jhwpe5lObd7GkuH1lPAFEcrnpbp+7v63SX5PsaeygVa4x1RGVjIxk2LIkZODrG
GL1gRfdQIfMZ4Eeb/TF8HmkYDNmp9CRZ2mb/EAKhn6kBMob/Os0dfTZlkUr3tSRlVZzfG5MdjaJg
FV5jVTj+8HICke/5Ohk3DA94scUssN2AqmOw8kODn/p4uDHPr9R9WTdTlFc0Knv3Jr+6Mc7NVTAY
wq+fI2MggR8xnCbG86XieAm1r1INABzngWmoqQ/MjxmN8C4LYrlWK3ZY4UwKbpI6KlyLOpfvI6Ov
mvv8JCqRitMzdqbpHCjF9gW142wdq/JFHgdevNSrwW/+xHAQ
"""
        self.par.update(
            {
                "PROJET": "Miss_Laplace",
                "DREF": 5,
                "FICHIER_SOL_INCI": "./Miss_Laplace.sol.inci",
                "FREQ_IMAG": 64.97465,
                "INST_FIN": 2.0,
                "LIST_FREQ": (1.023383,),
                "NB_MODE": 6,
                "OFFSET_MAX": 40,
                "OFFSET_NB": 400,
                "PAS_INST": 0.05,
                "PRECISION": 1e-10,
                "RFIC": 0.0,
                "SPEC_MAX": 0.075,
                "SPEC_NB": 16384,
                "SURF": "OUI",
                "Z0": -5.8,
            }
        )
        gen = MissCmdeGen(self.par, self.struct, self.fname)
        txt = gen.build()
        if self._debug:
            with open("/tmp/test05_miss03c_loop.in", "w") as f:
                f.write(txt)
        diff = self._diffcompress(refe, txt)
        assert diff.strip() == "", diff

    def test06_zzzz108c(self):
        """ISS using TABLE_CONTROL"""
        refe = """
eJzNU9tu1DAQffdXzDtqSFIV9aUPxh5nLeILvqxWfUGrbYCVSku3q4rPZ5xs6e4WBKJCIkosz8zJ
mePL6dBigEt6mvpcMMkTZ0mngMwHJ7NO4HkAwXvxwegYmZm7HqrXjz9UXx5ur1kXXPYMAFqgcjbI
lLYlPt2P6euYmHHjWa8D7rOsPn9l2ibsAk/aWQgoEryBc0hBc9v1CE1LL1MB32e0AkEilH7VWV0e
PKkb4CVzVtVT5hVlPI/HKBalM1zb8fcGRuW7Ocll74QpOuMBrH2CnRRYy2IKfLv+uB4mMA2S7fMy
XCQMGnN4ttTh27bgcR/fUhBbA5d1ieqnJVAcc1DPSO5p07UVmhRL33NiiXNoqp8QsLqimGo/gLM/
BPrf4XCBAoqKaSacTcHRUWVNy1FhuIPkcsIIkc8RzFxJSiQqQ6JPGy9BuSBevhNj++hVmGa5Sx1k
WcN410ojvFrerIax2zBhut695f1h57EgtVKb5Wo7QC5iEw1/LansyOb2eo8KrOvQMbRSbYY7FnVn
eQ/HBHVT3a8/3Ryq8y4mFj3RWnVRMnS/leGLi6NLT5dLzI7YViSpWjIRybs9dtQegQuBPU5+i6Ox
dsc1WXRnsJ2rTmH0UkMqyKbxF0Up+8IRdw56gdaHQ61zTboi/ncyrw5lSiTnCDRIvf+FVBp8Gb4D
RYFWCQ==
"""
        self.par.update(
            {
                "PROJET": "ZZZZ108C",
                "FREQ_MIN": 0.25,
                "FREQ_MAX": 50.0,
                "FREQ_PAS": 0.25,
                "OFFSET_MAX": 60,
                "OFFSET_NB": 600,
                "SURF": "OUI",
                "ALGO": "DEPL",
                "Z0": 0.0,
                "NBM_TOT": 291,
                "NBM_DYN": 291,
                "NBM_STA": 291,
                "TOUT_CHAM": "OUI",
                "_calc_impe": False,
                "_calc_forc": False,
                "_hasPC": True,
                "_nbPC": 3,
            }
        )
        gen = MissCmdeGen(self.par, self.struct, self.fname)
        txt = gen.build()
        if self._debug:
            with open("/tmp/test06_zzzz108c.in", "w") as f:
                f.write(txt)
        diff = self._diffcompress(refe, txt)
        assert diff.strip() == "", diff

    def test07_fdlv112e(self):
        """use ISSF + control points"""
        refe = """
eJzdV1mP4kYQfu9fUc+zGsKRifKyDx677bHia7vbaJWXyGMasGRs1kc0+fepbiAYMOCdzEpREIfp
46vjq66uegBLvrYNLCTkIJsWf/2srskDPPa9yAPOBOUGVrKQVfatlbizhmWWrjNZ1eC7nF/b24Vx
aEAZ2JY3n0ymVKPu9WhhI4sWLEMYvUBkN4NzImsq2VGcPFwuFq5glEQstGJXQGQwMA3P/GOnJmJ4
Mm1ahYJikyzPk1UfjELy56EHo58OKo82f5b5Xu9lVmRNVhbaFauqbLfyqgPRdBbGEQGAJ0DI2KfE
dgP1/5fuf/w4ffCbcnEVnJgvhh8Rz2W0q2m63mqkKKmSjUSf1bCArGjkqkoU7g2+iBsI6jBDuGEA
jJoClfwVBHONwPEoTKb4vtQS8gS2ypXqeVlJjJIiveGSM5k2o19iGpgULIpugZ9HY/2in8ZjMNTI
U3ckMrgam3TG+hxXl239uCiR5eKuKhf725PtUDdVq+NmqEncCn3DDbQ5E9ARsH8GmCrqyW+mrzjn
d0WrqLsj9ETc9CjuUYmbKY8qHxIumNFky0z+GMGzjmC08XFGLNuLXeSUhWeEzcCkHuYDPfp0HH2m
QvM97vLNY2bjyFFpOyvupQ212OpmGZ003mTa3jgAer25TqqVRGydmt4RAKRLPaFfBWUujdnFKZVv
jdKS3pI6gIITeVfOJmqfKOLT5Mj0AOSjAD714ffxxUmEC6Pq/hyZrpPNtsYclGYLtHFAZsBEZLqA
IRR5BtrG5xgoPSqQ8QhDA+f+WfgycGF0b52mJcnTNn+HCcg7NUHZcA8HkvYNtmWmHpGttCyaqsyH
Jpo+qWYYCBZito7dG7G1zNts8Z3hPNN4z2Wb5hicbQW5KgUG5XuEsZn8BiKMBeXAjTkFf25bOCBQ
URD4cf3IAjtk5v/zTHSCYFkWqRKrKXcqKW/dyTuBmloe2ewSq0LfY0UE2WYrF8mwm3eHFzvCgdga
gy4kFAFUA2gW5A+Jno/UvE9fJtHf7YHTbVW+5qj9vgCqkvROAXQpyfHCZ8P7+IjkfxXNWtaymxUO
Z0rrukyG11AHZS3XtpWREmJ1sAR+/VfC/5a9H5D/jrlPbT+aD0Ho0PBYN2jrX9+TwrqyAsvGPSeo
d8qL05ZnW9bNI1KQNTqcblXkPWVatiqSa2QR7jqB4cE5DePJSO37+DiOTk35t/cZiUIuCI+QzMD+
DLvC0vaNr5/Pb2rbNV/OLEzRhlFCTI6tm0ex58Ti0zB1ralbGq57jf0ltOue9j3HZNdoYLGs2osJ
aoGdEL8yaVmewuCHkhR/om5lOpTds/i504zrNvFvcCbg8g==
"""
        self.par.update(
            {
                "PROJET": "FDLV112E",
                "FREQ_MIN": 4.0,
                "FREQ_MAX": 5.0,
                "FREQ_PAS": 1.0,
                "Z0": 5.0,
                "SURF": "NON",
                "ALGO": "REGU",
                "OFFSET_MAX": 1000,
                "OFFSET_NB": 5000,
                "ISSF": "OUI",
                "MATER_FLUIDE": {
                    "RHO": 1000.0,
                    "CELE": 1500,
                    "AMOR_BETA": 0.0,
                    "DEMI_ESPACE": "OUI",
                },
                "NBM_DYN": 38,
                "NBM_STA": 6,
                "NBM_TOT": 44,
                "TOUT_CHAM": "NON",
                "_calc_impe": False,
                "_calc_forc": False,
                "_hasPC": True,
                "_nbPC": 3,
                "_hasSL": True,
            }
        )
        gen = MissCmdeGen(self.par, self.struct, self.fname)
        txt = gen.build()
        if self._debug:
            with open("/tmp/test07_fdlv112e.in", "w") as f:
                f.write(txt)
        diff = self._diffcompress(refe, txt)
        assert diff.strip() == "", diff

    def test08_fdlv113a(self):
        """idem fdlv112e + matériau homogène"""
        refe = """
eJy9Vtty2jAQfddX7HM6UAyh6atjy8QzxnZlmelbxjECPOML8aXTX+p39Me6kiEx4WLSaapBINba
o7O72l3fgCmemhqWAlIQdYO/86SqyA0MTg1yg0/cIoO1yEWZPDcCNStYJfEmEWUFczsIzul2YWbU
pQws01lo2kRXqDseDWQib8DUuX4SiLRP8BlP6lJ0iJOb483c5owSn3lmaHPwdQaG7hiPLU3EcERc
NxIFj42SNI3Wp2Ak0nzhOTD8vKc8zH4U6Y73KsmTOily5Yp1WTRbcdaBaDrzQp9Ytjs7pZ4Vy7PK
xHjQ5z5xbEa7TOLNViH5URllAn1SwRKSvBbrMpK4F+JBbJfTGdO57bnAqMHhC3wFzmzdnTkUtDF+
jllCGsFWukquV6XAW5DHF0x+c6bF6LeQugYFkwIObThSg34ajUCXkvFQ20s0DFrQyl52aaccVxVN
NVgWGMW8l8qRfnOgjn/SXmNIYHpz3XZbE0BFVa4HGn5NcN6Suc4pMK9lP92zn8DCl6LpcHp3dze5
RdEYFsrGyfBWm36ZjpTonnLljFHHPfLaBP+c/bjDfoxzQiwntM099058JmBQBxNXSV8sOs0VgpBZ
KHnlbCV5X37LzWa3HKjs/ini5sJNVvuNTVSuBWKrGvIeX2BJ6UTyAtgqbRJZZt6Ddyp9MFpljMkT
5RWk18MqaNs1bCQZeCEz6JHLe/4TaIc2bK2M0rhJVfbEmyjbVlg1YqSS1z0JROh3aoCkonDuiyZG
O6qmRHOqKysCOsli4hm4F3IaQKAvKMwXlokCDqENHKc9902wPLT0A0L8avyqwHBU2EQgybYCw3JV
MWudEM74TPKkSk2RFcAs20DHT+X8kPv0X8kzgS5u9qVmWxZPKZqy6zFlFPf0mONjZ453rzuvJUE1
lKe/uURdZNe0UOcAtadyHL52bIuqHtRllNQqUJe6puq2h9tVLDACv39dGQLiewHHcmc8dHs5Nu/m
UQYSgMjbb6huuUsS7P9+KNdqxdsVFkz/LKHd5aiSKpPvale36PPcJCJpU7LLzDSdlo7M3ANau6J/
rXffxK/nhVQeQ/4A5Lekkg==
"""
        self.par.update(
            {
                "PROJET": "FDLV113A",
                "FREQ_MIN": 1.0,
                "FREQ_MAX": 21.0,
                "FREQ_PAS": 20.0,
                "SURF": "NON",
                "RFIC": 0.5,
                "ISSF": "OUI",
                "MATER_SOL": {"E": 7.0e8, "RHO": 2500.0, "NU": 0.2, "AMOR_HYST": 0.0},
                "MATER_FLUIDE": {
                    "RHO": 1000.0,
                    "CELE": 150,
                    "AMOR_BETA": 0.0,
                    "DEMI_ESPACE": "NON",
                },
                "SOURCE_FLUIDE": {"POINT": (0.0, 0.0, 0.0)},
                "_hasSL": True,
            }
        )
        gen = MissCmdeGen(self.par, self.struct, self.fname)
        txt = gen.build()
        if self._debug:
            with open("/tmp/test08_fdlv113a.in", "w") as f:
                f.write(txt)
        diff = self._diffcompress(refe, txt)
        assert diff.strip() == "", diff

    def _diffcompress(self, refe, new):
        """Compare files without comment"""
        sref = remove_comments((test_utils.uncompress64(refe.encode())).decode())
        snew = remove_comments(new)
        diff = test_utils.difftxt(sref, snew)
        if diff and self._debug:
            with open("/tmp/miss_diff.refe", "w") as fref:
                fref.write(sref)
            with open("/tmp/miss_diff.new", "w") as fnew:
                fnew.write(snew)
        return diff


class TestMissInterface(unittest.TestCase):

    """test interface functions to create miss datafiles"""

    faster = "ZZZZ108B.aster"

    # unittest.skipIf(not osp.isfile(faster),   # decorator requires python 2.7
    # "requires %s" % faster)
    def test01_ext(self):
        """test creation of the .ext file"""
        if not osp.isfile(self.faster):
            return
        rdr = ResuAsterReader(nbgrp=2)
        data = rdr.read(self.faster)
        assert data.noeud_nb == 100, data.noeud_nb
        assert data.maille_nb_tot == 99, data.maille_nb_tot
        assert data.maille_nb == [96, 3], data.maille_nb
        assert data.mode_dyna_nb == 0, data.mode_dyna_nb
        assert data.mode_stat_nb == 291, data.mode_stat_nb


class TestMissCsolReader(unittest.TestCase):

    """test the reader of csol files"""

    fcsol = "ZZZZ108B.01.csol.a"

    # unittest.skipIf(not osp.isfile(faster),   # decorator requires python 2.7
    # "requires %s" % faster)
    def test01_ext(self):
        """test creation of the .ext file"""
        if not osp.isfile(self.fcsol):
            return
        reader = MissCsolReader(3, 201)
        lfreq, values = reader.read(self.fcsol)
        self.tab = tab = Table()
        tab["FREQ"] = lfreq
        for ipc, respc in enumerate(values):
            for iddl, comp in enumerate(("X", "Y", "Z")):
                lab = "PC_{}_{}_REEL".format(ipc + 1, comp)
                tab[lab] = respc.comp[iddl][0]
                lab = "PC_{}_{}_IMAG".format(ipc + 1, comp)
                tab[lab] = respc.comp[iddl][1]


if __name__ == "__main__":
    unittest.main()
