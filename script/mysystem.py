# Hello emacs, this is -*- python -*-
# $Id: mysystem.py,v 1.4 2001/08/03 14:05:38 andre Exp $
# André Rabello <Andre.Rabello@ufrj.br>

# This module controls how the program will interact with the system, creating
# directories and redirecting output.

# What will we need to proceed
import os
import sys
import os.path
import getpass
import getopt
import string
import time
import socket
import shutil

# Private library
import sdb #for configuration
import re #for perl regexp
import builder #for network building
import tester #for network testing
import analyzer #for network analysis

# This class defines the program configuration, machine parameters, etc.
class _Properties:
    def __init__(self):
        "Creates the run configuration"
        # Who am I talking to
        self.__owner=getpass.getuser()
        # The program name
        self.__program_name=os.path.basename(sys.argv[0])
        # Where I started
        self.__start_dir=os.getcwd()
        
    def start_dir(self):
        return self.__start_dir
    
    def home(self):
        "Go back to where we started"
	print self.prompt(), 'Changing to home dir (' + self.__start_dir + ')'
        os.chdir(self.__start_dir)
        
    def prompt(self):
        "Print a simple prompt that indicates the program name and user"
        # This will define what to print when called
        if len(self.__program_name) != 0:
            return '[%s:%s]' % (self.__owner, self.__program_name)
        else:
            return '[%s:interactive]' % self.__owner

    def error(self):
        "What to print in case of error"
        return self.prompt()+' *ERROR*'
        
# This class defines the methods for interacting with the system, given the
# configuration above
class System:
    "Defines all system operations"
    __config = None
    __wd = None
    
    def __init__(self,arguments=""):
        "Loads the run configuration"
        self.props=_Properties() # get current system properties
        if len(arguments) != 0:
            print self.props.prompt(), 'Processing input options...'
        else:
            print self.props.prompt(), 'Using default options...'

        # try to process input options
        try:
            self.__opts, self.__args = \
                         getopt.getopt(arguments,'c:b:hw:')
            #print self.__args
            if len(self.__args) != 0:
                raise InvalidArgument(self.__args[0])
        except getopt.error, error:
            print self.props.error(), 'getopt:', error
            print self.help()
            sys.exit(2)
        except InvalidArgument, error:
            print self.props.error(), 'invalid argument:', error
            print self.help()
            sys.exit(2)
        except:
            print "[MySystem] Unexpected error:", sys.exc_info()[0]
            raise

        # loads the default configuration
        self.__config = sdb.Parser('./default.sdb')

        # loads the configuration first if needed
        newopt = []
        for o, a in self.__opts:
            if o == '-c':
                if os.path.exists(a):
                    print self.props.prompt(),'Configuration overrided by',a
                    self.__config = sdb.Parser(a)
                else:
                    print self.props.error(),'config file: can\'t find',a
                    sys.exit(2)
                    
            else: #append to new set of options
                newopt.append((o,a))
                
        self.__opts = newopt
        
        # Start processing the options
        for o, a in self.__opts:
            if o == '-h':
                print self.help()
                sys.exit(2)
                
            elif o == '-b': #change where to put the results
                try:
                    self.__config_base_wd(a)
                except:
                    print "[MySystem] Unexpected error:", sys.exc_info()[0]

            elif o == '-w': #define an already existing working directory
                if os.path.exists(a+'/run.dat'): #if the directory is valid
                    self.__config_wd(a)
                else:
                    mesg='The directory you gave me is _not_ from a known run!'
                    mesg2="\n Or, it doesn't exist, which is worse..."
                    print self.props.error(), mesg+mesg2
                    sys.exit(2)

    def __config_wd(self,dir):
        "[Private] Configure base and working directory from another run"
        if dir[len(dir)-1] == '/': #eliminate such char
            dir = dir[:len(dir)-1]
        split_path = string.split(dir,'/')
        bwd = split_path[0]
        if len(split_path) < 2: #I need, at least, the base and working dirs
            mesg='The directory you gave me exists, but I cannot work'
            mesg2=' with the working directory only!'
            print self.props.error(), mesg+mesg2
            sys.exit(2)

        for i in range(1,len(split_path)-1):
            bwd = bwd + '/' + split_path[i]
        self.__config_base_wd(bwd)

        #Our last resource go to the working directory variable, we are
        # almost all set now...
        self.__wd = split_path[len(split_path)-1]
        
        #Configure the run time of the test
        stime = re.search("(?P<time>\d{7,20})",self.__wd)
        self.__starttime = time.localtime(int(stime.group('time')))

        #Give a warning message
        print self.props.prompt(),\
              'Base and Working directory configured successfuly!'

    def __config_base_wd(self,dir):
        "[Private] Configure the working directory"

        # First, check if the directory exists
        if os.path.exists(dir): #do nothing
            self.__config.set('/config/base_wd',dir)
            return

        if string.find(dir,'/') == 0:
            checkpath = ''
        else:
            checkpath = '.'
        for i in string.split(dir,'/'): #all but last directory
            checkpath = checkpath+ '/' + i
            if not os.path.exists(checkpath):
                os.mkdir(checkpath)
                print self.props.prompt(),'Created directory', checkpath

        #Now I can configure into DB
        match = re.match(os.environ['HOME']+'(?P<rest>.*)',dir)
        if match != None:
            current = match.group('rest')
            self.__config.set('/config/base_wd','@HOME@'+dir)
        else:
            current = os.getcwd()
            self.__config.set('/config/base_wd',dir)

        print self.props.prompt(),'Redirected output to base directory',dir

    def help(self):
        "Print a help message"
        print self.props.prompt(), 'Configuration Options'
        print self.props.prompt(), \
              '-c <string>     Configuration File (in python)'
        print self.props.prompt(), \
              '-b <string>     The base path for the output directory'
        print self.props.prompt(), \
              '-w <string>     The path for the output directory (leave empty!)'
        print self.props.prompt(), \
              '-h            Print this help message' 

    def chbwd(self):
        "Change directory to base working directory"
        base_wd = self.__config.get('/config/base_wd')
        if not os.path.exists(base_wd): #checks the base working diretory
            self.__config_base_wd(base_wd)
	print self.props.prompt(),'Changing to', base_wd, '...'
        os.chdir(base_wd)

    def chwd(self):
        "Change directory to working directory"

        # The following process will create a unique working
        # directory, based on the machine basename and the current time
        # in UNIX ticks. This can only be overriden _if_ and only if the
        # machine where you are running the programs have it's current
        # time modified (changed to past times). In such case, there's a
        # possibility of overriding the output directory, if you choose to run
        # these routines _exactly_ at the same instant as previous times and
        # with the same base working diretory. This shall be rare enough for
        # us _not_ to bother.
        if self.__wd == None: #we didn't set the working directory
            time.sleep(1) #for us to have different working directories
            self.__starttime = time.localtime(time.time()) #get current time
            host = string.split(socket.gethostname(),'.')[0] #the base hostname
            self.__wd = host+'@'+time.strftime('%s',self.__starttime)

        self.chbwd() #changes to base working directory
        bwd = self.__config.get('/config/base_wd')
	print self.props.prompt(),'Changing to', bwd+'/'+self.__wd, '...'
        if not os.path.exists(self.__wd): #checks the working directory
            os.mkdir(self.__wd)
            #copy the current configuration into the working directory
            self.__config.write(self.__wd+'/default.sdb')
            if os.path.exists('last'):
                os.remove('last')
            os.symlink(self.__wd,'last') #create a symbolic link from 'last'
        os.chdir(self.__wd)

    def chhd(self):
        "Change to home directory (config wrapper)"
        self.props.home()

    def build(self):
        "Will build the neural network based on configuration read"
        args = self.__config.get('/config/train')
        args.update(self.__config.get('/config/network'))
        b = builder.Builder(args,\
                            self.__config.get('/config/program/builder'),\
                            self.__config.get('/config/script/dryer'))
        self.chwd()
        b.build()
        print self.props.prompt(), "Network built"
        self.chhd()

    def run_test(self,section):
        "Defines how to test a set with the build network"
        args = self.__config.get(section) #set special details for this run
        args.update(self.__config.get('/config/run')) #add run time params.
        args.update(self.__config.get('/config/network')) #network params.
        t = tester.Tester(self.__config.get('/config/program/tester'),args)
        self.chwd()
        t.test()
        print self.props.prompt(), "Network tested for", section
        self.chhd()

    def test_validation_set(self):
        if self.__config.exists('/config/validation_set'):
            self.run_test('/config/validation_set')

    def test_test_set(self):
        self.run_test('/config/test_set')

    def test_train_set(self):
        self.run_test('/config/train_set')

    def analyze(self):
        print self.props.prompt(), "Analyzing results..."
        self.chwd()
        a = analyzer.Analyzer(self.__config.get("/analysis"))
        a.analyze()
        self.chhd()
        print self.props.prompt(), "Analysis finished."

    def force_compress(self):
        """Will force working directory compression.

        Even if configuration file says the contrary.
        """
        if self.__config.exists('/config/compress_wd'):
            self.__config.set('/config/compress_wd','yes')
        
    def compress_wd(self):
        "Will compact the working directory in a tar.gzipped file"
        #if it's told to do so...
        if self.__config.exists('/config/compress_wd') and \
           self.__config.get('/config/compress_wd') != 'yes':
            return

        else:
            #else, goes on compressing the working directory
            wd = self.__wd
            self.chbwd() #go to base working directory
            
            #get the time in ticks in which the test started
            volname = time.strftime( \
                    '"Net built at %A, %B %d of %Y at %H:%m:%S"', \
                    self.__starttime)
            volarg = '--label='+volname #set archive (volume) name or label
            runlist = '/bin/tar cfz '+wd+'.tar.gz '+wd+' '+volarg
            os.system(runlist)

            shutil.rmtree(wd) #remove the whole working directory
            print self.props.prompt(), \
                  "Working directory compressed and removed"
            self.chhd()

##############################################################################
# EXCEPTIONS                                                                 #
##############################################################################

# Defines the error to throw in case of the existance of workdir
class WorkDirExists:
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value + ' already exists!'
 
# Defines the error in case the user supplies extra arguments
class InvalidArgument:
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return 'Extra argument not allowed: ' + self.value





