#! /usr/bin/env python

import gcat_workflow.core.gcat_conf as gc
import gcat_workflow.core.run_conf as rc

def main(args):

    ###
    # set run_conf
    run_conf = rc.Run_conf(
        sample_conf_file = args.sample_conf_file,
        project_root = args.project_root,
        gcat_conf_file = args.gcat_conf_file,
    )
    #run_conf.analysis_type = args.analysis_type
    #run_conf.drmaa = False if args.disable_drmaa else True
    run_conf.runner =  args.runner
    run_conf.retry_count = args.retry_count

    # disable params
    # - args.multiprocess
    # - args.param_check
    
    ###
    # set gcat_conf and task parameter config data
    gcat_conf = gc.gcat_conf(conf = run_conf.gcat_conf_file)
    gcat_conf.software_version_set()
    
    if args.analysis_type == "dna":
        import gcat_workflow.dna.sample_conf as sc
        import gcat_workflow.dna.configure as configure

    elif args.analysis_type == "rna":
        import gcat_workflow.rna.sample_conf as sc
        import gcat_workflow.rna.configure as configure
    
    elif args.analysis_type == "germ":
        import gcat_workflow.germ.sample_conf as sc
        import gcat_workflow.germ.configure as configure
        
    sample_conf = sc.Sample_conf(run_conf.sample_conf_file)
    configure.main(gcat_conf = gcat_conf, run_conf = run_conf, sample_conf = sample_conf)

