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

# person_in_charge: mathieu.courtois at edf.fr

"""
Definition of the Timer object.
"""

import os
import time
from contextlib import contextmanager


from .i18n import translate as _


def _dtimes():
    """Return a dict of cpu, system and total times."""
    l_t = os.times()
    return {"cpu": (l_t[0], l_t[2]), "sys": (l_t[1], l_t[3]), "tot": l_t[4]}


def _conv_hms(t):
    """Convert a number of seconds in hours, minutes, seconds."""
    h = int(t / 3600)
    m = int(t % 3600) // 60
    s = (t % 3600) % 60
    return h, m, s


class Timer:
    """This class provides methods to easily measure time spent during
    different steps.
    Methods :
       Start : start a timer in mode 'INIT' ([re]start from 0) or 'CONT'
          (restart from last value).
       Stop  : stop a timer
    Attributes :
       timers : dict {
          timer_id : {
             'name'   : timer legend (=timer_id by default),
             'state'  : state,
             'cpu_t0' : initial cpu time,
             'cpu_dt' : spent cpu time,
             'sys_t0' : initial system time,
             'sys_dt' : spent system time,
             'tot_t0' : start time,
             'tot_dt' : total spent time,
             'num'    : timer number (to print timers in order of creation),
             'hide'   : boolean,
          },
          ...
       }
          state is one of 'start', 'stop'
    """

    MaxNumTimer = 9999999
    MaxSize = 500

    def __init__(self, add_total=True, format="as_run", maxlabel=None, limit=None, title=None):
        """Constructor"""
        # ----- initialisation
        self.timers = {}
        self.add_total = add_total
        if limit:
            self.MaxSize = limit
        self._oversize = False
        self._oversize_name = maxlabel or _("after %d timers")
        try:
            self._oversize_name = self._oversize_name % self.MaxSize
        except TypeError:
            pass
        if not format in ("as_run", "aster"):
            format = "as_run"

        if format == "as_run":
            self.fmtlig = (
                "   %(name)-33s  %(cpu_dt)9.2f  %(sys_dt)9.2f  %(cpu_sys)9.2f  %(tot_dt)9.2f"
            )
            self.fmtstr = "   %(title)-33s  %(cpu)9s  %(sys)9s  %(cpu+sys)9s  %(elapsed)9s"
            self.sepa = " " + "-" * 81
            self.TotalKey = _("Total time")
            self.d_labels = {
                "title": "",
                "cpu": _("cpu"),
                "sys": _("system"),
                "cpu+sys": _("cpu+sys"),
                "elapsed": _("elapsed"),
            }
        elif format == "aster":
            self.fmtlig = " * %(name)-24s : %(cpu_dt)10.2f : %(sys_dt)10.2f : %(cpu_sys)10.2f : %(tot_dt)10.2f *"
            self.fmtstr = (
                " * %(title)-24s : %(cpu)10s : %(sys)10s : %(cpu+sys)10s : %(elapsed)10s *"
            )
            self.sepa = " " + "*" * 80
            self.TotalKey = "TOTAL_JOB"
            self.d_labels = {
                "title": "COMMAND",
                "cpu": "USER",
                "sys": "SYSTEM",
                "cpu+sys": "USER+SYS",
                "elapsed": "ELAPSED",
            }
        if title:
            self.TotalKey = title

        self.total_key = id(self)
        if self.add_total:
            self.Start(self.total_key, name=self.TotalKey, num=self.MaxNumTimer)

    def Start(self, timer, mode="CONT", num=None, hide=None, name=None):
        """Start a new timer or restart one"""
        if len(self.timers) > self.MaxSize and name != self._oversize_name:
            if not self._oversize:
                self._oversize = True
                self.Start(self._oversize_name, num=num, name=self._oversize_name)
        name = name or str(timer)
        isnew = self.timers.get(timer) is None
        if not num:
            num = len(self.timers)
        if mode == "INIT":
            num = self.timers[timer]["num"]
        dico = _dtimes()
        if isnew or mode == "INIT":
            self.timers[timer] = {
                "name": name,
                "state": "start",
                "cpu_t0": dico["cpu"],
                "cpu_dt": 0.0,
                "sys_t0": dico["sys"],
                "sys_dt": 0.0,
                "tot_t0": dico["tot"],
                "tot_dt": 0.0,
                "num": num,
                "hide": hide,
            }
        elif mode == "CONT" and self.timers[timer]["state"] == "stop":
            self.timers[timer].update(
                {
                    "state": "start",
                    "cpu_t0": dico["cpu"],
                    "sys_t0": dico["sys"],
                    "tot_t0": dico["tot"],
                }
            )

    def Add(self, timer, cpu_dt=0.0, sys_dt=0.0, to_total=False):
        """Add dt values (hidden to os.times, for example under mpiexec) to a timer."""
        if len(self.timers) > self.MaxSize:
            return
        if self.timers.get(timer) is not None:
            self.timers[timer]["cpu_dt"] += cpu_dt
            self.timers[timer]["sys_dt"] += sys_dt
        if to_total and timer != self.total_key:
            self.Add(self.total_key, cpu_dt, sys_dt)

    def Stop(self, timer, hide=None):
        """Stop a timer"""
        if self.timers.get(timer) is None:
            if len(self.timers) > self.MaxSize:
                return
            self.timers[timer] = {
                "name": str(timer),
                "hide": hide,
                "state": "stop",
                "cpu_t0": 0.0,
                "cpu_dt": 0.0,
                "sys_t0": 0.0,
                "sys_dt": 0.0,
                "tot_t0": 0.0,
                "tot_dt": 0.0,
                "num": len(self.timers),
            }
        elif self.timers[timer]["state"] == "start":
            dico = _dtimes()
            self.timers[timer]["state"] = "stop"
            for i in range(len(dico["cpu"])):
                self.timers[timer]["cpu_dt"] += dico["cpu"][i] - self.timers[timer]["cpu_t0"][i]
            self.timers[timer]["cpu_t0"] = dico["cpu"]
            for i in range(len(dico["sys"])):
                self.timers[timer]["sys_dt"] += dico["sys"][i] - self.timers[timer]["sys_t0"][i]
            self.timers[timer]["sys_t0"] = dico["sys"]
            self.timers[timer]["tot_dt"] = (
                self.timers[timer]["tot_dt"] + dico["tot"] - self.timers[timer]["tot_t0"]
            )
            self.timers[timer]["tot_t0"] = dico["tot"]
            if hide is not None:
                self.timers[timer]["hide"] = hide

    def StopAndGet(self, timer, *args, **kwargs):
        """Stop a timer and return "delta" values."""
        self.Stop(timer, *args, **kwargs)
        if self.timers.get(timer) is None:
            return 0.0, 0.0, 0.0
        cpu_dt = self.timers[timer]["cpu_dt"]
        sys_dt = self.timers[timer]["sys_dt"]
        tot_dt = self.timers[timer]["tot_dt"]
        if len(self.timers) > self.MaxSize:
            del self.timers[timer]
        return cpu_dt, sys_dt, tot_dt

    def StopAndGetTotal(self):
        """Stop the timer and return total "delta" values."""
        return self.StopAndGet(self.total_key)

    def getsortedtimers(self):
        """Return timers list sorted by timer number."""
        lnum = [
            [timer["num"], timer]
            for timer in list(self.timers.values())
            if timer["hide"] is not True
        ]
        lnum = sorted(lnum, key=lambda x: x[0])
        return lnum

    def StopAll(self):
        """Stop all timers"""
        lk = list(self.timers.keys())
        if self.add_total:
            lk.remove(self.total_key)
        for timer in lk:
            self.Stop(timer)

    def __repr__(self):
        """Pretty print content of the timer.
        NB : call automatically StopAll
        """
        self.StopAll()
        if self.add_total:
            self.Stop(self.total_key)

        labels = self.fmtstr % self.d_labels
        out = [""]
        # get timers list and sort by 'num'
        lnum = self.getsortedtimers()
        if lnum:
            out.append(self.sepa)
            if self.add_total and labels:
                out.append(labels)
                out.append(self.sepa)
        for num, timer in lnum:
            d_info = timer.copy()
            d_info["cpu_sys"] = d_info["cpu_dt"] + d_info["sys_dt"]
            if self.add_total and num == self.MaxNumTimer and len(lnum) > 1:
                out.append(self.sepa)
            out.append(self.fmtlig % d_info)
        if lnum:
            out.append(self.sepa)
        out.append("")
        return os.linesep.join(out)

    @contextmanager
    def __call__(self, timer, mode="CONT", num=None, hide=None, name=None):
        """Context manager used to measure the time spent in a block."""
        self.Start(timer, mode, num, hide, name)
        yield
        self.Stop(timer, hide)


if __name__ == "__main__":
    chrono = Timer(format="aster")
    with chrono("Global"):
        chrono.Start("Compilation")
        chrono.Start("CALC_FONCTION")
        chrono.Start(23, name="CALC_FONCTION")
        time.sleep(0.4)
        chrono.Stop("Compilation")
        chrono.Stop(23)
        chrono.Start("Child")
    with chrono("Other"):
        time.sleep(0.15)
    # first CALC_FONCTION will be stop by 'repr()'
    print(chrono)
