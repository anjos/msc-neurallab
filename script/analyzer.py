# Hello emacs, this is -*- python -*-
# $Id: analyzer.py,v 1.1 2001/07/13 15:48:02 andre Exp $
# André Rabello <Andre.Rabello@ufrj.br>

# This module controls how the program analyze the results of training,
# testing and validating the neural network. The core of analysis is done
# by Matlab, therefore this shall only be a wrapper.

import tempfile
import os
import shutil
import string

_LOCATION = os.environ['HOME'] + '/work/neural/tools/analyze.m'

class Analyzer:
    "This class is responsible for analyzing the outputs of a particular run"

    def __init__(self,opts):
        "Initializes the analyzer"
        self.__flags = opts
        self.__matsh = opts['program']['matsh']
        self.__analyzer = opts['program']['analyzer']
        if not os.path.exists(self.__analyzer):
            print '[ANALYZER] Cannot locate analyzer =>', self.__analyzer
            sys.exit(0)

    def change(self, what):
        "Changes the parameters required by the second argument"
        for k in what.keys():
            if self.__flags.has_key(k):
                self.__flags[k] = what[k]
            else:
                raise ParameterDoesNotExist(k)

    def parse_efficiencies(self):
        "Parses the instance efficiencies, transforming them into matlab code"

        effs = self.__flags['efficiency']
        refs = []

        for i in effs.keys(): #get the references
            current = {'line':'','label':'','file':''}
            for j in effs[i].keys(): #get all features
                
                if j == 'line_type':
                    lt = effs[i][j]
                    if lt == 'filled':
                        current['line'] = current['line'] + '-';
                    elif lt == 'dotted':
                        current['line'] = current['line'] + '.';
                    elif lt == 'dashdot':
                        current['line'] = current['line'] + '.-';
                    elif lt == 'dashed':
                        current['line'] = current['line'] + '--';
                    else:
                        raise OptionNotImplemented(lt,j)

                elif j == 'color':
                    c = effs[i][j]
                    if c == 'red':
                        current['line'] = 'r' + current['line']
                    elif c == 'blue':
                        current['line'] = 'b' + current['line']
                    elif c == 'green':
                        current['line'] = 'g' + current['line']
                    elif c == 'black':
                        current['line'] = 'k' + current['line']
                    elif c == 'magenta':
                        current['line'] = 'm' + current['line']
                    elif c == 'yellow':
                        current['line'] = 'y' + current['line']
                    elif c == 'white':
                        current['line'] = 'w' + current['line']
                    else:
                        raise OptionNotImplemented(c,j)

                elif j == 'label' or j == 'file':
                    current[j] = effs[i][j]

                else:
                    raise ParameterDoesNotExist(j)

            refs.append(current)
            
        del current

        #From refs, create the code as a structure containing the data
        init = []
        for i in refs: #build matlab structures with those
            current = "struct('name', '"+i['file']+\
                      "','label', '"+i['label']+\
                      "','line', '"+i['line']+"')"
            init.append(current)

        line = ""
        for i in init: #build the matlab commmand line
            if len(line) == 0:
                line = 'efflist = ['+i+','
            else:
                line = line+i+','
        #when we end, we have an extra comma ',' at the end, erase it!
        line = line[0:len(line)-1]+'];\n'
        # Finally, we return the code itself...
        return line

    def parse_filenames(self):
        "Parse filenames that are necessary for data processing"

        files = self.__flags['file']
        line = []

        for i in files.keys(): #get the references
            line.append(i+" = '"+files[i]+"';\n")
        return line

    def analyze(self):
        "Will run the analysis itself, parsing arguments and calling Matlab"
        init_code = [self.parse_efficiencies()]
        init_code.extend(self.parse_filenames())

        #Build a temporary that will be called to initilize Matlab...
        print '[ANALYZER] Initializing Matlab from temporary script'
        tfilename = tempfile.mktemp()
        tfile = open(tfilename,'w')
        tfile.writelines(init_code)
        tfile.close()
        initfile = self.__flags['init']
        shutil.copyfile(tfilename,initfile)
        os.remove(tfilename) #delete, as we no longer need it...
        initscript = string.replace(initfile,'.m','')
        script = string.replace(os.path.basename(self.__analyzer),'.m','')
        command = self.__matsh+' "'+initscript+';'+\
                  'addpath '+os.path.dirname(self.__analyzer)+";"+\
                  'eval(\''+script+'\');"'
        #Now I can call matlab with the required parameters
        print command
        os.system(command)

#############################################################################
# Exceptions
#############################################################################
class ParameterDoesNotExist:
    def __init__(self, key):
        self.key = key
    def __str__(self):
        return "[ANALYZER] Parameter "+self.key+" does not exist on config!"

class OptionNotImplemented:
    def __init__(self,opt,what):
        self.opt = opt
        self.what = what
    def __str__(self):
        return '[ANALYZER] Option "' + self.opt + \
               '" not implemented when scanning "'+ self.what + '"'

