"""
Root package for all executable CSB client programs.

Follow these simple steps to write a new CSB app:

    1. Create the app module in the C{csb.apps} package.
    
    2. Create a main class and derive it from L{csb.apps.Application}. You need
    to implement the L{csb.apps.Application.main()} abstract method - this is
    the app's entry point. You have the L{csb.apps.Application.args} object at
    your disposal.
    
    3. Create MyAppRunner class, derived from csb.apps.AppRunner. You need to
    implement the following methods and properties:
    
        - property L{csb.apps.AppRunner.target} -- just return YourApp's class
        - method L{csb.apps.AppRunner.command_line()} -- make an instance of
        L{csb.apps.ArgHandler}, define your command line parameters on that
        instance and return it
        - optionally, override L{csb.apps.AppRunner.initapp(args)} if you need
        to customize the instantiation of the main app class, or to perform
        additional checks on the parsed application C{args} and eventually call
        C{YourApp.exit()}. Return an instance of your app at the end 
    
    4. Make it executable::
        if __name__ == '__main__':
            MyAppRunner().run()
    
See L{csb.apps.helloworld} for a sample implementation.    
"""

import os
import re
import sys
import argparse
import traceback

from abc import ABCMeta, abstractmethod, abstractproperty


class ExitCodes(object):
    """
    Exit code constants.
    """
    CLEAN = 0
    USAGE_ERROR = 1
    CRASH = 99
    
class AppExit(Exception):
    """
    Used to signal an immediate application exit condition (e.g. a fatal error),
    that propagates down to the client, instead of forcing the interpreter to
    close via C{sys.exit()}.
    
    @param message: exit message
    @type message: str
    @param code: exit code (see L{ExitCodes} for common constants)
    @type code: int
    @param usage: ask the app runner to print also the app's usage line
    @type usage: bool
    """
    
    def __init__(self, message='', code=0, usage=False):
        
        self.message = message
        self.code = code
        self.usage = usage

class Application(object):
    """
    Base CSB application class.
    
    @param args: an object containing the application arguments
    @type args: argparse.Namespace
    """
    __metaclass__ = ABCMeta
    
    USAGE = ''
    HELP = ''
    
    def __init__(self, args, log=sys.stdout):
        
        self.__args = None
        self.__log = log
                
        self.args = args
        
    @property
    def args(self):
        """
        The object containing application's arguments, as returned by the
        command line parser.
        """
        return self.__args
    @args.setter
    def args(self, args):
        self.__args = args
    
    @abstractmethod
    def main(self):
        """
        The main application hook.
        """
        pass
    
    def log(self, message, ending='\n'):
        """
        Write C{message} to the logging stream and flush it.
        
        @param message: message
        @type message: str        
        """
        
        self.__log.write(message)
        self.__log.write(ending)
        self.__log.flush()
        
    @staticmethod
    def exit(message, code=0, usage=False):
        """
        Notify the app runner about an application exit.

        @param message: exit message
        @type message: str
        @param code: exit code (see L{ExitCodes} for common constants)
        @type code: int
        @param usage: advice the client to show the usage line
        @type usage: bool        
        
        @note: you re not supposed to use C{sys.exit()} for the same purpose.
               It is L{AppRunner}'s responsibility to handle the real system
               exit, if the application has been started as an executable.
               Think about your app being executed by some Python client as a
               regular Python class, imported from a module -- in that case you
               only want to ask the client to terminate the app, not to kill
               the whole interpreter.   
        """
        raise AppExit(message, code, usage)        
            
class AppRunner(object):
    """
    A base abstract class for all application runners. Concrete sub-classes
    must define their corresponding L{Application} using the L{self.target}
    property and must customize the L{Application}'s command line parser using
    L{self.command_line()}.
    
    @param argv: the list of command line arguments passed to the program. By
                 default this is C{sys.argv}.
    @type argv: tuple of str 
    """
    __metaclass__ = ABCMeta
        
    def __init__(self, argv=sys.argv):
        
        self.module = argv[0]
        self.program = os.path.basename(self.module)
        self.args = argv[1:]
    
    @abstractproperty
    def target(self):
        """
        Reference to the concrete {Application} class to run. This is
        an abstract property that couples the current C{AppRunner} to its
        corresponding L{Application}.

        @rtype: type (class reference)
        """
        return Application 
    
    @abstractmethod
    def command_line(self):
        """
        Command line factory: build a command line parser suitable for the
        application.
        This is a hook method that each concrete AppRunner must implement.
        
        @return: a command line parser object which knows how to handle
        C{sys.argv} in the context of the concrete application. See the
        documentation of L{ArgHandler} for more info on how to define command
        line arguments.
        
        @rtype: L{ArgHandler}
        """
        # null implementation (no cmd arguments):
        return ArgHandler(self.program)
    
    def initapp(self, args):
        """
        Hook method that controls the instantiation of the main app class.
        If the application has a custom constructor, you can adjust the
        app initialization by overriding this method.
        
        @param args: an object containing the application arguments
        @type args: argparse.Namespace
        
        @return: the application instance
        @rtype: L{Application}
        """
        app = self.target
        return app(args)

    def run(self):
        """
        Get the L{self.command_line()} and run L{self.target}. Ensure clean
        system exit. 
        """
        try:
            app = self.target
            cmd = self.command_line()
                    
            try:           
                assert issubclass(app, Application)
                assert isinstance(cmd, ArgHandler)
                            
                args = cmd.parse(self.args)
                app.USAGE = cmd.usage
                app.HELP = cmd.help
    
                self.initapp(args).main()
                
            except AppExit as ae:
                if ae.usage:
                    AppRunner.exit(ae.message, code=ae.code, usage=cmd.usage)
                else:
                    AppRunner.exit(ae.message, code=ae.code)
    
            except SystemExit as se:                            # this should never happen, but just in case 
                AppRunner.exit(se.message, code=se.code)
                        
        except Exception:
            message = '{0} has crashed. Details: \n{1}'.format(self.program, traceback.format_exc())
            AppRunner.exit(message, code=ExitCodes.CRASH)
        
        AppRunner.exit(code=ExitCodes.CLEAN)            
            
    @staticmethod
    def exit(message='', code=0, usage='', ending='\n'):
        """
        Perform system exit. If the exit C{code} is 0, print all messages to
        STDOUT, else write to STDERR.
        
        @param message: message to print
        @type message: str
        @param code: application exit code
        @type code: int   
        """
        
        ending = str(ending or '')
        message = str(message or '')
        stream = sys.stdout
        
        if code > 0:
            message = 'E#{0} {1}'.format(code, message)
            stream = sys.stderr
        
        if usage:
            stream.write(usage.rstrip(ending))            
            stream.write(ending)
        if message:
            stream.write(message)            
            stream.write(ending)
        
        sys.exit(code)

class ArgHandler(object):
    """
    Command line argument handler.
    
    @param program: (file)name of the program, usually sys.argv[0]
    @type program: str
    @param description: long description of the application, shown in help
                        pages. The usage line and the parameter lists are
                        generated automatically, so no need to put them here.
    @type description: str
    
    @note: a help argument (-h) is provided automatically. 
    """
    
    SHORT_PREFIX = '-'
    LONG_PREFIX = '--'
    
    class Type(object):
        
        POSITIONAL = 1
        NAMED = 2
    
    def __init__(self, program, description=''):
        
        self._argformat = re.compile('^[a-z][a-z0-9_-]*$', re.IGNORECASE)
        self._optformat = re.compile('^[a-z0-9]$', re.IGNORECASE)
        
        self._program = program
        self._description = description
        
        self._parser = argparse.ArgumentParser(prog=program, description=description)
        
    def _add(self, kind, name, shortname, *a, **k):
        
        args = []
        kargs = dict(k)
                    
        if shortname is not None:
            if not re.match(self._optformat, shortname):
                raise ValueError('Invalid short option name: {0}.'.format(shortname))

            if kind == ArgHandler.Type.POSITIONAL:
                args.append(shortname)
            else:                     
                args.append(ArgHandler.SHORT_PREFIX + shortname)

        if name is not None or kind == ArgHandler.Type.POSITIONAL:
            if not re.match(self._argformat, name):
                raise ValueError('Malformed argument name: {0}.'.format(name))
            
            if kind == ArgHandler.Type.POSITIONAL:
                args.append(name)
            else:
                args.append(ArgHandler.LONG_PREFIX + name)

        assert len(args) in (1, 2)   
        args.extend(a)                        
        
        self.parser.add_argument(*args, **kargs)        
        
    def add_positional_argument(self, name, type, help, choices=None):
        """
        Define a mandatory positional argument (an argument without a dash).
        
        @param name: name of the argument (used in help only)
        @type name: str
        @param type: argument data type
        @type type: type (type factory callable)
        @param help: help text
        @type help: str
        @param choices: list of allowed argument values
        @type choices: tuple
        """
        self._add(ArgHandler.Type.POSITIONAL, name, None,
                  type=type, help=help, choices=choices)

    def add_array_argument(self, name, type, help, choices=None):
        """
        Same as L{self.add_positional_argument()}, but allow unlimited number
        of values to be specified on the command line.
        
        @param name: name of the argument (used in help only)
        @type name: str
        @param type: argument data type
        @type type: type (type factory callable)
        @param help: help text
        @type help: str
        @param choices: list of allowed argument values
        @type choices: tuple
        """
        self._add(ArgHandler.Type.POSITIONAL, name, None,
                  type=type, help=help, choices=choices, nargs=argparse.ONE_OR_MORE)        

    def add_boolean_option(self, name, shortname, help, default=False):
        """
        Define an optional switch (a dashed argument with no value).
        
        @param name: long name of the option (or None)
        @type name: str, None
        @param shortname: short (single character) name of the option (or None)
        @type shortname:str, None
        @param help: help text
        @type help: str
        @param default: default value, assigned when the option is omitted. 
                        If the option is specified on the command line, the
                        inverse value is assigned  
        @type default: bool       
        """
        if not help:
            help = ''
        help = '{0} (default={1})'.format(help, default)
        
        if default:
            action = 'store_false'
        else:
            action = 'store_true'
                     
        self._add(ArgHandler.Type.NAMED, name, shortname,
                  help=help, action=action, default=bool(default))
        
    def add_scalar_option(self, name, shortname, type, help, default=None, choices=None, required=False):
        """
        Define a scalar option (a dashed argument that accepts a single value).
        
        @param name: long name of the option (or None)
        @type name: str, None
        @param shortname: short (single character) name of the option (or None)
        @type shortname: str, None
        @param type: argument data type
        @type type: type (type factory callable)        
        @param help: help text
        @type help: str
        @param default: default value, assigned when the option is omitted
        @param choices: list of allowed argument values
        @type choices: tuple
        @param required: make this option a named mandatory argument
        @type required: bool      
        """
        if not help:
            help = ''
        if default is not None:
            help = '{0} (default={1})'.format(help, default)             
         
        self._add(ArgHandler.Type.NAMED, name, shortname,
                  type=type, help=help, default=default, choices=choices, required=required)        

    def add_array_option(self, name, shortname, type, help, default=None, choices=None, required=False):
        """
        Define an array option (a dashed argument that may receive one
        or multiple values on the command line, separated with spaces).

        @param name: long name of the option (or None)
        @type name: str, None
        @param shortname: short (single character) name of the option (or None)
        @type shortname: str, None
        @param type: argument data type
        @type type: type (type factory callable)        
        @param help: help text
        @type help: str
        @param choices: list of allowed argument values
        @type choices: tuple
        @param required: make this option a named mandatory argument
        @type required: bool                   
        """
        if not help:
            help = ''
        if default is not None:
            help = '{0} (default={1})'.format(help, default)           
         
        self._add(ArgHandler.Type.NAMED, name, shortname,
                  nargs=argparse.ZERO_OR_MORE, type=type, help=help, choices=choices, required=required)
        
    def parse(self, args):
        """
        Parse the command line arguments.
        
        @param args: the list of user-provided command line arguments --
                     normally sys.argv[1:]
        @type args: tuple of str        
        
        @return: an object initialized with the parsed arguments
        @rtype: argparse.Namespace
        """
        try:
            return self.parser.parse_args(args)
        except SystemExit as se:
            if se.code > 0:
                raise AppExit('Bad command line', ExitCodes.USAGE_ERROR)
            else:
                raise AppExit(code=ExitCodes.CLEAN)                
    
    @property
    def parser(self):
        return self._parser
    
    @property
    def usage(self):
        return self.parser.format_usage()

    @property
    def help(self):
        return self.parser.format_help()    
