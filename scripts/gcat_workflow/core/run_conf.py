#! /usr/bin/env python

import os
class Run_conf(object):
    """
    class for job related parameters
    """

    def __init__(self, sample_conf_file = None, 
                       project_root = None, 
                       gcat_conf_file = None):

        self.sample_conf_file = sample_conf_file
        self.gcat_conf_file = gcat_conf_file 
        self.project_root = ""
        os.makedirs(project_root, exist_ok = True)
        if project_root.startswith("/"):
            self.project_root = project_root
        else:
            import subprocess
            import tempfile
            temp_sh = tempfile.mktemp()
            temp_out = tempfile.mktemp()
            try:
                f = open(temp_sh, "w")
                f.write("cd %s\n" % (project_root))
                f.write("pwd > %s" % (temp_out))
                f.close()
                subprocess.call(["bash", temp_sh])
                self.project_root = open(temp_out).read().rstrip("\n")
                #print(self.project_root)
            finally:
                os.remove(temp_sh)
                os.remove(temp_out)
                
