"""
Patch for MSVC manifest embedding to handle intermittent file locking issues.

This module patches the ccroot.link_task.exec_mf() method to add retry logic
when MT.exe fails with "Access is denied" errors, which is a common issue on
Windows due to antivirus scanners or other processes temporarily locking
the executable file.
"""

import time
from waflib import Logs, Utils
from waflib.Tools import ccroot


# Store the original exec_mf method
_original_exec_mf = ccroot.link_task.exec_mf


def exec_mf_with_retry(self):
    """
    Execute MT.exe to embed manifest with retry logic for file locking issues.

    This patched version will retry up to 3 times with exponential backoff
    if MT.exe fails with an access denied error (exit code 31).
    """
    if not self.env.MT:
        return 0

    manifest = None
    for out_node in self.outputs:
        if out_node.name.endswith('.manifest'):
            manifest = out_node.abspath()
            break
    else:
        return 0

    mode = ''
    for x in Utils.to_list(self.generator.features):
        if x in ('cprogram', 'cxxprogram', 'fcprogram', 'fcprogram_test'):
            mode = 1
        elif x in ('cshlib', 'cxxshlib', 'fcshlib'):
            mode = 2

    Logs.debug('msvc: embedding manifest in mode %r', mode)

    lst = [] + self.env.MT
    lst.extend(Utils.to_list(self.env.MTFLAGS))
    lst.extend(['-manifest', manifest])
    lst.append('-outputresource:%s;%s' % (self.outputs[0].abspath(), mode))

    # Retry logic for file locking issues
    max_retries = 3
    retry_delay = 0.5  # seconds
    ret = 1  # Default to error if somehow no attempts are made

    for attempt in range(max_retries):
        ret = super(ccroot.link_task, self).exec_command(lst)

        if ret == 0:
            return 0

        # Exit code 31 typically indicates MT.exe failed (often due to file locking)
        if ret == 31 and attempt < max_retries - 1:
            delay = retry_delay * (2 ** attempt)  # exponential backoff
            Logs.warn('msvc: MT.exe failed (attempt %d/%d), retrying in %.1f seconds...' %
                     (attempt + 1, max_retries, delay))
            time.sleep(delay)
            continue

        # If we've exhausted retries or it's a different error, return the error code
        if attempt == max_retries - 1:
            Logs.error('msvc: MT.exe failed after %d attempts' % max_retries)

        return ret

    return ret


def configure(conf):
    """
    Configure function to patch the manifest embedding with retry logic.
    This is called automatically when the module is loaded.
    """
    # Only apply the patch on Windows with MSVC/Intel compilers
    if Utils.is_win32:
        Logs.info('Applying MSVC manifest retry patch for file locking issues')
        ccroot.link_task.exec_mf = exec_mf_with_retry

