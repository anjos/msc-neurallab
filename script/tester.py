# Hello emacs, this is -*- python -*-
# $Id: tester.py,v 1.2 2001/08/03 14:07:28 andre Exp $
# André Rabello <Andre.Rabello@ufrj.br>

# This module defines the tester class, which is responsible for
# testing the new neural network built. As parameters, you should pass the
# correspondent dictionary of terms (containing configuration)

import sys
import os
import re
import string

class _Converter:
    "Defines a nice way to load the net configuration from JETNET dumps"

    def __init__(self,file='net.dat'):
        "Will read in the file containg the network definitions in JETNET"
        fp = open(file,'r')
        self.__save = fp.readlines()
        fp.close()

    def convert(self,file=0):
        "Will create a raw network configuration file that is used by testers"
        
        i = 0;
        layer = [] #temp holder for layer configuration
        config = {}
        while 1:
            if i == len(self.__save):
                break #we are finished!
            
            l = self.__save[i] #get next line
            i = i + 1 #prepare to get the other one

            # test 1: network configuration
            reob = re.match('\s*(?P<int>[\d]{1,3})\snodes.*\s(?P<ln>\d).*',l)
            if reob != None: #got net properties
                layer.insert(0,int(reob.group('int')))
                node = ''
                if layer[0] == 1:
                    node = 'node'
                else:
                    node = 'nodes'
                print '[Converter] Layer',int(reob.group('ln')),\
                      'has',layer[0],node
                continue #stop

            # test 2: input layer configuration
            reob = re.match('\s*.*layer\s(?P<start>\d).*and\s(?P<end>\d).*',l)
            if reob != None: # got weights
                lname = reob.group('start')+'->'+reob.group('end')
                print '[Converter] Processing data between',lname
                config[lname] = {'in': layer[int(reob.group('start'))],\
                                  'out': layer[int(reob.group('end'))],\
                                  'weights': [],\
                                  'bias': []}

                i = i + 1 #a blank line
                for j in range(0,config[lname]['in']):
                    weights = ''
                    while len(string.strip(self.__save[i])) != 0:
                        weights = weights + ' ' + self.__save[i]
                        i = i + 1
                    weights = string.strip(weights)
                    config[lname]['weights'].append(weights)
                    i = i + 1 #a blank line between valid ones

                #convert all char lines into floats
                config[lname]['weights'] = self.__convert_matrix(\
                    config[lname]['weights'])

                #go to the threshold field
                i = i + 2
                bias = ''
                while i < len(self.__save) and \
                              len(string.strip(self.__save[i])) != 0:
                    bias = bias + ' ' + self.__save[i]
                    i = i + 1
                bias = string.strip(bias)
                config[lname]['bias'].append(bias)

                #convert all char lines into floats
                config[lname]['bias'] = self.__convert_matrix(\
                    config[lname]['bias'])
                
                continue #stop

        #print config
        accline = [] #accumulator of lines to write
        accline.append(str(config['0->1']['in']) + ' ' + \
                       str(config['0->1']['out']) + '\n')
        #write-out for layers 0->1 and 1->2 only (1 hidden layer only)
        for i in range(0,len(config['0->1']['weights'][0])): #for each column
            l = '' #empty string
            for j in range(0,len(config['0->1']['weights'])): #and each line
                l = l + config['0->1']['weights'][j][i] + ' '
            l = l + config['0->1']['bias'][0][i] + '\n'
            accline.append(l)
        for i in range(0,len(config['1->2']['weights'][0])): #for each column
            l = '' #empty string
            for j in range(0,len(config['1->2']['weights'])): #and each line
                l = l + config['1->2']['weights'][j][i] + ' '
            l = l + config['1->2']['bias'][0][i] + '\n'
            accline.append(l)

        #now I can do the file stuff
        if file != 0:
            fp = open(file,'w')
        else:
            fp = sys.stdout
        fp.writelines(accline)
        if file != 0:
            fp.close()

    def __convert_matrix(self, lst):
        "gets a list of lines with N values per line and put in matrix form"
        tmp = []
        for i in range(0,len(lst)):
            tmp.append(string.split(lst[i]))
            for j in range(0,len(tmp[i])): #convert into floats
                tmp[i][j] = string.strip(tmp[i][j])
                #tmp[i][j] = float(tmp[i][j]) #convertion if needed
        return tmp

class Tester:
    "Defines the neural network tester"

    def __init__(self, program, opts):
        self.__tester = program
        if not os.path.exists(self.__tester):
            print '[TESTER] Cannot locate tester =>', self.__tester
            #sys.exit(0)

        #opts is a dictionary will all data in it
        self.__flags = opts

    def change(self, what):
        "Changes the parameters required by the second argument"
        for k in what.keys():
            if self.__flags.has_key(k):
                self.__flags[k] = what[k]
            else:
                raise ParameterDoesNotExist(k)

    def test(self):
        "Executes the neural testing with the parameters indicated"
        args = []
        for k in self.__flags.keys(): #build arguments
            if k == 'file':
                args.append('--input-file='+self.__flags[k])
                if not os.path.exists(self.__flags[k]):
                    raise CannotFindFile(self.__flags[k])
            elif k == 'activation':
                if self.__flags[k] == 'builtin':
                    args.append('--act-builtin')
            elif k == 'hint':
                args.append('--hint='+self.__flags[k])
            elif k == 'relevance':
                if self.__flags[k] == 'true':
                    args.append('--relevance')
            elif k == 'importance':
                if self.__flags[k] == 'true':
                    args.append('--importance')
            elif k == 'output':
                args.append('--output-file='+self.__flags[k])
            elif k == 'rawnet':
                args.append('--net-file='+self.__flags[k])
                self.__rawNetFile = self.__flags[k]
            elif k == 'net':
                self.__netFile = self.__flags[k]
            else:
                print '[TESTER] Option "'+k+'" not used by tester!'

        # First we have to convert the network file produced by JETNET
        # into something more usable
        if not os.path.exists(self.__rawNetFile):
            print '[TESTER] Translating network genearated by Builder...'
            c = _Converter(self.__netFile) #read
            c.convert(self.__rawNetFile) #write
            print '[TESTER] Network translated!'
            
        #Converts list of arguments into a single string
        runlist = self.__tester
        for i in args:
            runlist = runlist + ' ' + i
        os.system(runlist)
            
#############################################################################
# Exceptions
#############################################################################
class ParameterDoesNotExist:
    def __init__(self, key):
        self.key = key
    def __str__(self):
        return "[TESTER] Parameter "+self.key+" does not exist on config!"

class CannotFindFile:
    def __init__(self,file):
        self.file = file
    def __str__(self):
        return "[TESTER] Can't find file "+self.file
