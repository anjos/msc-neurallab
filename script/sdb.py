# Hello emacs, this is -*- python -*-
# $Id: sdb.py,v 1.4 2001/08/02 01:25:13 andre Exp $
# André Rabello <Andre.Rabello@ufrj.br>

# This module defines a set of functions that can lookup a configuration file,
# returning its key values and section parameters. It's very generic and you
# can use it in your own applications. Unless you want deeper sections, you
# won't have to alter this source file.
# In order to use this, just create the sdb file like this
#> [section]
#> int var1 = value1
#> float var2 = value2
#> - subsection
#> string var3 = value3
#> -- subsubsection
#> float var4 = value4
#> float var5 = value5
#> - other_subsection
#> string var6 = value6
# For strings, you can reference other strings you've created in the file
# and use environment variables. The syntax is the following:
# @HOME@/x/y/z => is converted /your/home/directory/x/y/z
# %var3% => is converted to value3 (as defined above)
#
# End then import this module
#>>> import sdb
# And instantiate a configuration
#>>> config = sdb.Parser('config.sdb') #put the filename you choose
# And you can lookup a variable name
#>>> config.get('/section/var1')
#>>> config.get('/section/subsection/subsubsection/var5')
# And, of course, you can alter the values and re-dump the configuration
# into a file with sdb.Parser.write()
#>>> config.set('/section/var1',newvalue)
#>>> config.write('new-file.sdb')
# or into your screen
#>>> config.write()

import re # get regexp stuff (à la perl)
import string
import os
import sys
import types # get the types names StringType FloatType IntType NoneType

_prompt = '[SDB]'
_version = 0.1

class _Reader:
    "Defines the (configurable) Simple Database Reader"
    __dict = {}
    __lines = []
    __count = 0
    
    def __init__(self,file,dict):
        self.__dict = dict
        try:
            f = open(file, 'r')
            self.__lines = f.readlines()
            # print _prompt, file , 'contains' , len(self.__lines) , 'lines'
            f.close()
        except:
            print "[SDB] ERROR: Unexpected exception:", sys.exc_info()[0]
            raise

    def next(self):
        "Return the next tag"
        # self.__count = self.__count + 1
        if len(self.__lines): # prevent from popping zero length lists
            cline = self.__lines.pop(0)
        else:
            return None # it's the end of file

        cline = string.strip(cline) # get rid of excessive whitespaces
        # print '[',self.__count,']',cline
        rem = re.split('\s*\#',cline,1) # rem.pop(0) has data
        data = rem.pop(0) # get what I need
        for e in self.__dict.keys():
            x=re.match(self.__dict[e],data)
            if x:
                return (e,x)
        # if you got here you have nothing :)
        return ('ignore',None); # this line is nothing, ignore

class Parser:
    "This defines the class that will actually read the configuration file"
    __dict = None #my dictionary of valid terms
    __reader = None #my db reader
    __params = None #where I keep the parameters
    __temp = None #a temporary space
    __mytype = 'int|float|string' # the defined types
    __myname = '\w+[\w\-\.]*' # naming convention
    __mychars = '\w\-\.\$\:\%\/,\{\}@áéíúâêôçüà'
    __myvalue = '['+__mychars+']?[\s'+__mychars+']*'
    
    def __init__(self,file):
        "Loads configuration parameters"
        # Define our dictionary for database commands
        myname = self.__myname
        mychars = self.__mychars
        myvalue = self.__myvalue
        self.__dict = {'section': '^\s*\[(?P<name>'+myname+')\]', \
                       'subsection': '^\s*-\s?(?P<name>'+myname+')', \
                       'subsubsection': '^\s*--\s?(?P<name>'+myname+')', \
                       'var': '^\s*(?P<type>'+self.__mytype+')?\s+(?P<name>'+myname+')\s*\=\s*(?P<value>'+myvalue+')'}
        print '[SDB] Trying to load',file,'...'
        self.__reader = _Reader(file,self.__dict)
        
        try: #trying to parse the configuration file
            self.__params = self.__parse()
        except (DepthNotAllowed,VersionDoNotMatch), retval:
            print retval
            return None
        except:
            print '[SDB] Unexpected exception:', sys.exc_info()[0]
            raise
        print '[SDB] Configuration file read succesfuly!'
    
    def __replace_environ(self,value):
        #Eliminate @ from possible chars
        mychars = string.replace(self.__mychars,'\@','')
        name = re.findall('\@(?P<name>['+mychars+']+)\@',value)
        for i in name:
            try:
                rep = os.environ[i]
            except KeyError:
                print "[SDB] Can't find environment variable $"+i
                print "[SDB] Error when processing", value
                return None
            
            # print 'Replacing environment variable', i, 'by', rep
            value = string.replace(value,'@'+i+'@',rep)

        return value

    def __replace_variable(self,value):
        "Replaces defined variables and returns the final version of it"
        #Eliminate % from possible chars
        mychars = string.replace(self.__mychars,'\%','')
        name = re.findall('\%(?P<name>['+mychars+']+)\%',value)
        for i in name:
            try:
                rep = self.get(i)
            except KeyError:
                print "[SDB] Can't find defined variable %"+i+"%"
                print "[SDB] Error when processing", value
                return None
            
            # print 'Replacing defined variable', val, 'by', rep
            value = string.replace(value,'%'+i+'%',rep)

        return value

    def __replace_references(self,item):
        "Returns a version of the given dictionary, with all strings replaced\
         by the equivalent environment or defined variable. This procedure\
         only works with strings, so don't expect to reference floats and\
         integers."
        
        if type(item) == types.DictType: #continue with this
            for i in item.keys():
                if type(item[i]) == types.DictType: #recurse
                    item[i] = self.__replace_references(item[i])
                elif type(item[i]) == types.StringType: #apply changes
                    item[i] = self.__replace_environ(item[i])
                    item[i] = self.__replace_variable(item[i])
                    
        elif type(item) == types.StringType: #replace data
            item = self.__replace_environ(item)
            item = self.__replace_variable(item)
            
        return item
        
    def __check_version(self,version):
        "Checks if version matches"
        current = version
        if current != _version:
            raise VersionDoNotMatch(_version,current)

    def __parse(self):
        "Loads the configuration file using the class Reader"
        # print _prompt, 'Scanning...'
        p = {} # empty dictionary
        level = 0
        level_name = 'root'
        sublevel_name = 'root'
        subsublevel_name = 'root'
        while 1:
            token = self.__reader.next()
            if token == None: #file is over!
                break
            
            if token[0] == 'ignore':
                continue #comments or empty lines

            if token[0] == 'section':
                level_name = token[1].group('name')
                p[level_name] = {}
                level = 1 # I'm inside this section
                continue

            if token[0] == 'subsection':
                if level >= 1:
                    sublevel_name = token[1].group('name')
                    p[level_name][sublevel_name]= {}
                    level = 2 # I'm inside this subsection
                    continue                    
                else:
                    raise DepthNotAllowed('subsection',level)

            if token[0] == 'subsubsection':
                if level >= 2:
                    subsublevel_name = token[1].group('name')
                    p[level_name][sublevel_name][subsublevel_name]= {}
                    level = 3 # I'm inside this subsection
                    continue
                else:
                    raise DepthNotAllowed('subsubsection',level)

            if token[0] == 'var':
                # Will accept any level here but 0, includding root!
                varname = token[1].group('name')
                type = token[1].group('type')
                if type == None: # defaults to string
                    type = 'string'

                # Now we do the conversions...
                value = token[1].group('value')
                if type == 'int': # needs conversion
                    value = int(value)
                if type == 'float': # needs conversion
                    value = float(value)

                if level == 0:
                    raise DepthNotAllowed('var',level) # variables only in sections and under

                # now, what we _can_ set
                if level == 1:
                    p[level_name][varname] = value
                    continue
                if level == 2:
                    p[level_name][sublevel_name][varname] = value
                    continue
                if level == 3:
                    p[level_name][sublevel_name][subsublevel_name][varname] = value
                    continue

            # gets the next token
            token = self.__reader.next()

        self.__check_version(p['sdb']['version'])
        # print _prompt, "Configuration file was read successfuly!"
        return p

    def issection(self,var):
        "Defines whether a variable is a section or a database-variable"
        if type(var) == types.DictType:
            return 1
        else:
            return 0

    def __recursive_write(self,path):
        "This writes the variables to a string and repasses the rest"
        # print 'Analyzing path =>', path
        # print 'Current stream is ', self.__temp
        if self.issection(self.__rawget(path)):
            split_path = string.split(path,'/') 
            level = len(split_path)-1
            current = split_path.pop()
            if not len(current) == 0: # eliminates the case where path == '/'
                if level == 1: #section
                    self.__temp.append('\n[' + current + ']\n')
                if level == 2: #subsection
                    self.__temp.append('\n- ' + current + '\n')
                if level == 3: #subsubsection
                    self.__temp.append('\n-- ' + current + '\n')

            for dir in self.__rawget(path).keys():
                if path == '/': #fix
                    newpath = path+dir
                else:
                    newpath = path+'/'+dir

                if not self.issection(self.__rawget(newpath)): #write now!
                    value = self.__rawget(newpath)
                    varname = string.split(newpath,'/').pop()
                    vartype = str(type(value))
                    vartype = re.match("\<type\ \'(?P<type>"+self.__mytype+")'>",vartype).group('type')
                    value = str(value) #convert into string format
                    self.__temp.append(vartype+' '+varname+' = '+value+'\n')
                else:
                    self.__temp.append(self.__recursive_write(newpath))

    def write(self,name=0):
        "Writes the configuration to stdout or to a file, fast"
        self.__temp = []
        self.__recursive_write('/')
        # fix a strange behaviour on __recursive_write()
        def isNotNone(var):
            return None != var
        self.__temp = filter(isNotNone,self.__temp)
        if name != 0:
            fp = open(name,'w')
        else:
            fp = sys.stdout
        fp.writelines(self.__temp)
        fp.close
        self.__temp = None #clear
                                
    def exists(self,path):
        "Finds out whether some /path/you/want exists or not inside the tree."
        # Syntax is /the/path/you/want/varname
        split_path = string.split(path,'/')
        split_path.pop(0)
        temp=self.__params
        if len(split_path[0]) == 0: #in case we give just '/'
            return 1 #true
        for i in split_path: #all other cases
            try:
                temp=temp[i]
            except (KeyError):
                return 0 #false

        #In the case all tries didn't fail, we got a winner!
        return 1 #true

    def __rawget(self,path):
        """Locates the configuration data inside the configuration tree

        As in self.get(), but won't self.__replace_references!
        """
        # Syntax is /the/path/you/want/varname
        split_path = string.split(path,'/')
        split_path.pop(0)
        temp=self.__params
        if len(split_path[0]) == 0: #in case we give just '/'
            return temp
        for i in split_path: #all other cases
            try:
                temp=temp[i]
            except (KeyError):
                print 'The path you required does not exist =>',path
                raise #re-throw the exception

        return temp

    def get(self,path):
        "Locates the configuration data inside the configuration tree"
        # Syntax is /the/path/you/want/varname
        split_path = string.split(path,'/')
        split_path.pop(0)
        temp=self.__params
        if len(split_path[0]) == 0: #in case we give just '/'
            return self.__replace_references(temp)
        for i in split_path: #all other cases
            try:
                temp=temp[i]
            except (KeyError):
                print 'The path you required does not exist =>',path
                raise #re-throw the exception

        return self.__replace_references(temp)

    def set(self,path,value):
        "This will change the value of a variable, in the case it is not a section"
        if self.issection(path): #can't change that!
            raise CannotChangeSection(path)

        if type(value) != type(self.get(path)): #can't substitute by that!
            raise TypeChangeNotAllowed(type(self.get(path)),type(value))
        
        p = string.split(path,'/')
        p.pop(0) #not of interest
        if len(p) == 4: #subsubsection case
            self.__params[p[0]][p[1]][p[2]][p[3]] = value
        if len(p) == 3: #subsection case
            self.__params[p[0]][p[1]][p[2]] = value
        if len(p) == 2: #section case
            self.__params[p[0]][p[1]] = value
        
        #variable value altered
        return value

##############################################################################
# EXCEPTIONS                                                                 #
##############################################################################

# Defines the error to throw in case of no section
class DepthNotAllowed:
    def __init__(self,what,level):
        self.what = what
        self.level = level
    def __str__(self):
        return 'Cannot create ' + self.what + 'in level' + self.level    
 
# If the configuration file do not match SDB version
class VersionDoNotMatch:
    def __init__(self,capable,current):
        self.n = capable
        self.c = current
    def __str__(self):
        return 'Cannot read configuration with version different than', self.n, '... You have', self.c

# If you attempt to change a whole section
class CannotChangeSection:
    def __init__(self,path):
        self.p = path
    def __str__(self):
        return 'Cannot change a whole section! => [', self.p, ']'

# If you try to substitute type x by y
class TypeChangeNotAllowed:
    def __init__(self,typex,typey):
        self.x = str(typex)
        self.y = str(typey)
    def __str__(self):
        return 'Cannot change variable',self.x,'by',self.y
