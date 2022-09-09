#! /usr/bin/env python

import os
import datetime
import subprocess
import stat

class Runner(object):
    def __init__(self, singularity_script, qsub_option, log_dir, max_task, retry_count):
        self.qsub_option = qsub_option
        self.retry_count = retry_count
        self.singularity_script = os.path.abspath(singularity_script)
        #self.singularity_script = singularity_script
        self.jobname = os.path.basename(self.singularity_script).replace(".sh", "").replace("singularity_", "")
        #self.log_dir = os.path.abspath(log_dir)
        self.log_dir = log_dir
        self.max_task = max_task
        
    def task_exec(self):
        pass

class Drmaa_runner(Runner):
    def task_exec(self):
        import drmaa
    
        s = drmaa.Session()
        s.initialize()
         
        jt = s.createJobTemplate()
        jt.jobName = self.jobname
        jt.outputPath = ':' + self.log_dir
        jt.errorPath = ':' + self.log_dir
        jt.nativeSpecification = self.qsub_option
        jt.remoteCommand = self.singularity_script
        os.chmod(self.singularity_script, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP)

        returncode = 0
        returnflag = True
        if self.max_task == 0:
            for var in range(0, (self.retry_count+1)):
                jobid = s.runJob(jt)
                returncode = 0
                returnflag = True
                now = datetime.datetime.now()
                date = now.strftime("%Y-%m-%d %H:%M:%S")
                print ("Your job has been submitted with id: " + jobid + " at Date/Time: " + date)
                retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
                now = datetime.datetime.now()
                date = now.strftime("%Y-%m-%d %H:%M:%S")
                print ("Job: " + str(retval.jobId) + ' finished with status: ' + str(retval.hasExited) + ' and exit status: ' + str(retval.exitStatus) + " at Date/Time: " + date)
                returncode = retval.exitStatus
                returnflag = retval.hasExited
                if returncode == 0 and returnflag: break
            s.deleteJobTemplate(jt)
            s.exit()

        else:
            joblist = s.runBulkJobs(jt,1,self.max_task,1)
            all_jobids = []
            for var in range(0, (self.retry_count+1)):
                if len(all_jobids) > 0:
                    joblist = all_jobids
                    all_jobids = []
                returncode = 0
                returnflag = True
                now = datetime.datetime.now()
                date = now.strftime("%Y-%m-%d %H:%M:%S")
                print ('Your job has been submitted with id ' + str(joblist) + " at Date/Time: " + date)
                s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
                for curjob in joblist:
                    print ('Collecting job ' + curjob)
                    retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
                    now = datetime.datetime.now()
                    date = now.strftime("%Y-%m-%d %H:%M:%S")
                    print ("Job: " + str(retval.jobId) + ' finished with status: ' + str(retval.hasExited) + ' and exit status: ' + str(retval.exitStatus) + " at Date/Time: " + date)
                    
                    if retval.exitStatus != 0 or not retval.hasExited:
                        returncode = retval.exitStatus
                        returnflag = retval.hasExited
                        if var == self.retry_count: break
                        jobId_list = retval.jobId.split(".")
                        taskId = int(jobId_list[1])
                        all_jobids.extend(s.runBulkJobs(jt,taskId,taskId,1))
                   
                if returncode == 0 and returnflag: break
            s.deleteJobTemplate(jt)
            s.exit()

        if returncode != 0 or not returnflag: 
            raise RuntimeError("Job: " + str(retval.jobId)  + ' failed at Date/Time: ' + date)

class Qsub_runner(Runner):
    def task_exec(self):
        qsub_commands = ['qsub', '-sync', 'yes', '-N', self.jobname]
        if self.max_task != 0:
            qsub_commands.extend(['-t', '1-'+str(self.max_task)+':1'])

        qsub_options = []
        if type(self.qsub_option) == type(""):
            qsub_options += self.qsub_option.split(' ')
            if '' in qsub_options:
                qsub_options.remove('')
        if len(qsub_options) > 0:
            qsub_commands += qsub_options

        returncode = subprocess.call(qsub_commands + [
            "-o", self.log_dir,
            "-e", self.log_dir,
            self.singularity_script
        ])

        if returncode != 0: 
            raise RuntimeError("The batch job failed.")

# non-support array job
class Slurm_runner(Runner):
    def task_exec(self):
        qsub_commands = ['sbatch', '--wait']
        #if self.max_task != 0:
        #    qsub_commands.extend(['-t', '1-'+str(self.max_task)+':1'])

        qsub_options = []
        if type(self.qsub_option) == type(""):
            qsub_options += self.qsub_option.split(' ')
            if '' in qsub_options:
                qsub_options.remove('')
        if len(qsub_options) > 0:
            qsub_commands += qsub_options

        log_o_path = "%s/%s.o" % (self.log_dir, self.jobname) + "%j"
        log_e_path = "%s/%s.e" % (self.log_dir, self.jobname) + "%j"

        #print(qsub_commands + ["-o", log_o_path, "-e", log_e_path, "-J", self.jobname, self.singularity_script])
        returncode = subprocess.call(qsub_commands + [
            "-o", log_o_path, 
            "-e", log_e_path, 
            "-J", self.jobname,
            self.singularity_script
        ])

        if returncode != 0: 
            raise RuntimeError("The batch job failed.")

# non-support array job
class Bash_runner(Runner):
    def task_exec(self):
        bash_commands = ['bash']

        qsub_options = []
        #if type(self.qsub_option) == type(""):
        #    qsub_options += self.qsub_option.split(' ')
        #    if '' in qsub_options:
        #        qsub_options.remove('')

        if len(qsub_options) > 0:
            bash_commands = bash_commands + qsub_options

        log_file = self.log_dir + "/" + os.path.splitext(os.path.basename(self.singularity_script))[0] + ".log"
        returncode = -1
        with open(log_file, "w") as f:
            returncode = subprocess.call(bash_commands + [
                self.singularity_script
            ], stdout=f, stderr=f)
        if returncode != 0: 
            raise RuntimeError("The batch job failed.")

def main(args):
    import yaml
    conf = yaml.safe_load(open(args.conf))
    if args.interval > 0:
        import time
        import random
        time.sleep(int(args.interval*random.random())) 
    if conf["runner"] == "qsub":
        runner = Qsub_runner(args.script, conf["qsub_option"], conf["log_dir"], conf["max_task"], conf["retry_count"])
    elif conf["runner"] == "slurm":
        runner = Slurm_runner(args.script, conf["qsub_option"], conf["log_dir"], conf["max_task"], conf["retry_count"])
    elif conf["runner"] == "bash":
        runner = Bash_runner(args.script, conf["qsub_option"], conf["log_dir"], conf["max_task"], conf["retry_count"])
    else: 
        runner = Drmaa_runner(args.script, conf["qsub_option"], conf["log_dir"], conf["max_task"], conf["retry_count"])
    runner.task_exec()
