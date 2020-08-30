#!/usr/bin/env python
"""
This is a small script that does the equivalent of multicrab.
"""
import os
import sys
from optparse import OptionParser
import json
import copy
from multiprocessing import Process

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS --sampleType TYPE --era ERA --subEra SUBERA]"
             "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
             "\nUse multicrab -h for help")

    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = '',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-w', '--workArea',
                      dest = 'workArea',
                      default = 'crab',
                      help = "work area directory. Default: 'crab'.",
                      metavar = 'WAD')

    parser.add_option('-o', '--crabCmdOpts',
                      dest = 'crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD",
                      metavar = 'OPTS')

    parser.add_option('-e', '--era',
                      dest = 'era',
                      default = 'Run2018',
                      help = "Era to process: 'Run2018' (default), 'Run2017', 'Run2016'.",
                      metavar = 'ERA')

    parser.add_option('-u', '--subEra',
                      dest = 'subEra',
                      default = 'all',
                      help = "Sub-era to process: 'all' (default), custom (e.g. 'Run2016B').",
                      metavar = 'SUBERA')

    parser.add_option('-s', '--sampleType',
                      dest = 'sampleType',
                      default = 'all',
                      help = "Samples to process: 'all' (default), 'data', 'mc'. Only allowed if subEra='all'.",
                      metavar = 'TYPE')

    parser.add_option('-t', '--storageSite',
                      dest = 'storageSite',
                      default = 'CERN',
                      help = "Storage site: 'CERN' (default, note CERNBOX requires auth.), 'FNAL'.",
                      metavar = 'STORAGE')

    (options, arguments) = parser.parse_args()

    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options


def main():

    options = getOptions()

    doData = options.sampleType in ['all', 'data']
    doMC = options.sampleType in ['all', 'mc']

    era = options.era

    # The submit command needs special treatment.
    if options.crabCmd == 'submit':

        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------

        from CRABClient.UserUtilities import config, getUsername
        config = config()

        config.General.workArea = options.workArea
        config.General.transferOutputs = True
        config.General.transferLogs = True

        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = 'muonAnalyzer_cfg.py'
        #config.JobType.maxMemoryMB = 4000
        #config.JobType.numCores = 1

        config.Data.splitting = 'Automatic'

        config.Data.publication = False
        #config.Data.ignoreLocality = True

        if options.storageSite == 'FNAL':
            config.Site.storageSite = 'T3_US_FNALLPC'
            config.Data.outLFNDirBase = '/store/user/%s/TnP_ntuples/%s/' % (getUsername(), era)
            #config.Data.outLFNDirBase = '/store/group/lpcmetx/iDM/Muon_TnP/%s/' % era
        elif options.storageSite == 'CERN':
            # Note: CERNBOX write access from CRAB requires authorization
            # See https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3FAQ#Can_I_send_CRAB_output_to_CERNBO
            # and to request auth.: https://cern.service-now.com/service-portal?id=sc_cat_item&name=request-map-dn-to-gridmap&se=CERNBox-Service
            config.Site.storageSite = 'T2_CH_CERNBOX'
            config.Data.outLFNDirBase = '/store/user/%s/TnP_ntuples/%s/' % (getUsername(), era)

        #--------------------------------------------------------

        with open('../data/dataset_db.json', 'r') as db:
            data = json.load(db)

            samples = {}
            try:
                if options.subEra == 'all':
                    samples = data[era]
                else:
                    samples = dict({options.subEra: data[era][options.subEra]})
            except:
                print "Error! Requested era+sub-era is likely not valid. Please check argument."
                sys.exit()

        for sub_era, input_dataset in samples.items():
            
            isData = 'Run' in sub_era
            isRun2018D = '2018D' in sub_era

            if isData and not doData: continue
            if not isData and not doMC: continue

            if isData:
                #config.Data.splitting = 'LumiBased'
                #config.Data.unitsPerJob = 100
                if '2018' in era:
                    config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
                elif '2017' in era:
                    config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'
                elif '2016' in era:
                    config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt'

            config.JobType.pyCfgParams = ['isFullAOD={}'.format(True), 'isMC={}'.format(not isData), 'isRun2018D={}'.format(isRun2018D), 'numThreads=1', 'era={}'.format(era)]
            config.JobType.numCores = 1

            config.Data.inputDataset = input_dataset
            config.JobType.allowUndistributedCMSSW = True
            config.General.requestName = 'muonAnalyzer_' + era + '_' + sub_era
            #config.Data.outputDatasetTag = sample

            # If we need to pull input files from a list file instead of CRAB:
            # config.Data.userInputFiles = open(basedir + sample + '.list').readlines()

            # Submit.
            def submit(config, options):
                try:
                    print "Submitting for input dataset %s with options %s" % (input_dataset, options.crabCmdOpts)
                    crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
                except HTTPException as hte:
                    print "Submission for input dataset %s failed: %s" % (input_dataset, hte.headers)
                except ClientException as cle:
                    print "Submission for input dataset %s failed: %s" % (input_dataset, cle)

            # Need to submit using multiprocessing module because of CRAB issue with different configs
            p = Process(target=submit, args=(config,options,))
            p.start()
            p.join()

    # All other commands can be simply executed.
    elif options.workArea:

        for dir in os.listdir(options.workArea):
            projDir = os.path.join(options.workArea, dir)
            if not os.path.isdir(projDir):
                continue
            # Execute the crab command.
            msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
            print "-"*len(msg)
            print msg
            print "-"*len(msg)
            try:
                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers)
            except ClientException as cle:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle)

if __name__ == '__main__':
    main()
