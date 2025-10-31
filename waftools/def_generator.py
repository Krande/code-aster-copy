#!/usr/bin/env python3
"""
WAF tool for automatic DEF file generation during the build process.

This tool integrates .def file generation into the WAF build system,
ensuring that .def files are created before linking DLLs on Windows with MSVC.
"""

import os
import subprocess
from pathlib import Path
from waflib import Task, TaskGen, Utils
from waflib.Configure import conf


@conf
def find_dumpbin(conf):
    """Find dumpbin.exe in the system PATH."""
    conf.find_program('dumpbin', var='DUMPBIN', mandatory=False)


def configure(conf):
    """Configure the def_generator tool."""
    if Utils.is_win32:
        conf.find_dumpbin()


@TaskGen.feature('generate_def')
@TaskGen.before_method('process_source')
def generate_def_file(self):
    """
    Generate .def file from compiled object files before linking.

    This method is called during the build process to generate the .def file
    needed for exporting symbols from Windows DLLs.
    """
    # Get the library name and type
    lib_name = getattr(self, 'def_lib_name', self.name)
    lib_type = getattr(self, 'def_lib_type', 'fc')  # fc, cpp, or c

    # Determine the script to use
    script_map = {
        'fc': 'msvc/def_gen_fc.py',
        'cpp': 'msvc/def_gen_cpp.py',
        'c': 'msvc/def_gen_c.py',
    }

    script_path = self.bld.path.find_node(script_map.get(lib_type, script_map['fc']))
    if not script_path:
        self.bld.fatal(f'DEF generation script not found: {script_map.get(lib_type)}')

    # Output .def file path
    def_file = self.path.find_or_declare(f'msvc/{lib_name}.def')

    # Create the task for generating the .def file
    task = self.create_task('generate_def_task', src=script_path, tgt=def_file)
    task.lib_name = lib_name
    task.lib_type = lib_type
    task.env = self.env


class generate_def_task(Task.Task):
    """Task for generating .def files from object files."""

    color = 'BLUE'

    def run(self):
        """Execute the .def file generation."""
        # Determine build directory
        build_dir = self.generator.bld.variant_dir

        # Construct command to run the def generation script
        cmd = [
            self.env.PYTHON[0],
            self.inputs[0].abspath(),
            '--build-dir', build_dir,
            '--output', self.outputs[0].abspath(),
        ]

        # Run the command
        return self.exec_command(cmd)

    def keyword(self):
        return 'Generating'

    def __str__(self):
        return f'{self.outputs[0].nice_path()}'


@TaskGen.extension('.def')
def process_def_file(self, node):
    """Process .def files for linking."""
    # Add .def file to link flags
    self.link_task.env.append_value('LINKFLAGS', f'/DEF:{node.abspath()}')

