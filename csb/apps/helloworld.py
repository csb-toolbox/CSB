"""
This is a CSB HelloWorld dummy application.
"""

import sys
import csb.apps


class ExitCodes(csb.apps.ExitCodes):
    
    BAD_TEXT = 4
    
    
class AppRunner(csb.apps.AppRunner):
    
    @property
    def target(self):
        return HelloWorldApp
    
    def command_line(self):

        text =  "Hello World"
        
        cmd = csb.apps.ArgHandler(self.program, 'This program prints "Hello World".')
        
        cmd.add_scalar_option('text', 't', str, 'The text to print', default=text)
        cmd.add_boolean_option('upper', 'u', 'Print in upper case', default=False)
        
        return cmd
    
    def initapp(self, args):
        
        app = self.target
        if len(args.text) < 3:
            app.exit("Enter at least a few words", code=ExitCodes.BAD_TEXT, usage=True)
        
        return app(args)
    
    
class HelloWorldApp(csb.apps.Application):
    
    def main(self):
        
        if self.args.upper:
            text = self.args.text.upper()
        else:
            text = self.args.text
        
        self.log(text)
        self.log('HW: done.')        


def main():
    AppRunner(sys.argv).run()
    
    
if __name__ == '__main__':
    main()
