# Hello emacs, this is -*- python -*-
# $Id: builder.py,v 1.1 2001/07/13 15:48:02 andre Exp $
# André Rabello <Andre.Rabello@ufrj.br>

# This module defines the builder class, which is responsible for
# building the new neural network. As parameters, you should pass the
# correspondent dictionary of terms (containing configuration)

import sys
import os
import tempfile
import shutil

tempfile.tempdir = '/tmp' #sets the default temporary directory

class Builder:
    "Defines the neural network builder"

    def __init__(self, opts, program, dryer):
        self.__dryer = dryer
        self.__builder = program
        if not os.path.exists(self.__builder):
            print '[BUILDER] Cannot locate builder =>', self.__builder
            sys.exit(0)

        #opts is a dictionary will all data in it
        self.__flags = opts

    def change(self, what):
        "Changes the parameters required by the second argument"
        for k in what.keys():
            if self.__flags.has_key(k):
                self.__flags[k] = what[k]
            else:
                raise ParameterDoesNotExist(k)

    def build(self):
        "Executes the neural training with the parameters indicated"
        args = []
        for k in self.__flags.keys(): #build arguments
            if k == 'file':
                args.append('--train-file='+self.__flags[k])
                if not os.path.exists(self.__flags[k]):
                    raise CannotFindFile(self.__flags[k])
            elif k == 'test':
                args.append('--test-file='+self.__flags[k])
                if not os.path.exists(self.__flags[k]):
                    raise CannotFindFile(self.__flags[k])
            elif k == 'learning_rate':
                args.append('--learn-rate='+str(self.__flags[k]))
            elif k == 'momentum':
                args.append('--momentum='+str(self.__flags[k]))
            elif k == 'lr_decay':
                args.append('--lr-decay='+str(self.__flags[k]))
            elif k == 'batch':
                args.append('--batch='+str(self.__flags[k]))
            elif k == 'steps':
                args.append('--maxsteps='+str(self.__flags[k]))
            elif k == 'epoch':
                args.append('--epoch='+str(self.__flags[k]))
            elif k == 'runfile':
                args.append('--run-file='+self.__flags[k])
            elif k == 'config_file':
                args.append('--config-file='+self.__flags[k])
            elif k == 'effic_file':
                args.append('--eff-file='+self.__flags[k])
            elif k == 'input_dimension':
                args.append('--no-inputs='+str(self.__flags[k]))
            elif k == 'hidden_dimension':
                args.append('--no-hidden='+str(self.__flags[k]))
            elif k == 'network_saving':
                if self.__flags[k] == 'SP':
                    args.append('--use-sp')
                elif self.__flags[k] == 'Area':
                    args.append('--use-area')
                else:
                    raise InvalidParamValue(k,self.__flags[k])
            elif k == 'net': #this is for after training...
                self.__networkfile_name = self.__flags[k]
            else:
                print '[BUILDER] Config option','"'+k+'"','not used by Builder'

        #Converts list of arguments into a single string
        runlist = self.__builder
        for i in args:
            runlist = runlist + ' ' + i
        os.system(runlist)
        #now we have to fix the network file output by the builder
        self.dry('fort.10')
        print '[BUILDER] Renaming fort.10 as', self.__networkfile_name+'...'
        os.rename('fort.10',self.__networkfile_name)

    def dry(self,file):
        "Will readin the file, dry it, rewrite it and load the needed lines"
        tfile = tempfile.mktemp()
        print '[BUILDER] Working on temporary',tfile
        #gawk makes it easier to dry out the input file...
        os.system('/bin/gawk --file='+self.__dryer+' '+file+' > '+tfile)
        print '[BUILDER] Saving', tfile, 'as', file
        shutil.copyfile(tfile,file)
        os.remove(tfile)
        print '[BUILDER] File',file,'was dried out!'

            
#############################################################################
# Exceptions
#############################################################################
class ParameterDoesNotExist:
    def __init__(self, key):
        self.key = key
    def __str__(self):
        return "[BUILDER] Parameter "+self.key+" does not exist on config!"

class CannotFindFile:
    def __init__(self,file):
        self.file = file
    def __str__(self):
        return "[BUILDER] Can't find file "+self.file

class NetworkOutOfSpec:
    def __init__(self,nlay):
        self.nlay = nlay
    def __str__(self):
        return "[BUILDER] Can't work with more than "+self.nlay+"layers!"

class InvalidParamValue:
    def __init__(self,param,value):
        self.p = param
        self.v = value
    def __str__(self):
        return '[BUILDER] Parameter "'+self.p+'" _cannot_ take value "'+self.v+'"'
