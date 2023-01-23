"""
The prog ("program") module contains standard objects commonly used
in Python programs written by Lee. It was started by Lee on 080815.

"""


def get_help(name, description, options, usage, notes):
    """
    Creates and returns a standard help message that describes the program
    being run. The object needs to be passed four parameters: name, 
    description, usage, and notes, which are self-explanatory. Here is some 
    standard code to generate the parameters in the dependent program:

    name = '''
        programName.py
        '''
    description = '''
        First line of description...
        ...second line of description...
        '''
    options = '''
        -a --Description of option a.
        -b --Description of option b.
        '''
    usage = '''
        programName.py -a <whatisit?> -b <whatisit?>
        '''
    notes = '''
        First line of special things to note...
        '''

    Usage: help_message(name, description, options, usage, notes)
    """

    message = "\nNAME%s\nDESCRIPTION%s\nOPTIONS%s\nUSAGE%s\nNOTES%s" % (name, description, options, usage, notes)
    return message


def get_options(args, flags):
    """
    Checks and collects command line options for the program that calls it.
    Prints an error message if options are inappropriate or missing. Program
    must pass in expected options as a string, for example, 'hab:cdef:'.

    Usage: get_options(flags, message)
    """

    import sys
    import getopt

    try:
        options, xarguments = getopt.getopt(args, flags)
    except getopt.GetoptError:
        print '''
    Error: You used an unknown option or are missing an argument to an option. 
    Use the '-h' option for more information.
        '''
        sys.exit(0)
    return options

def check_help(options, message):
    """
    Checks command line arguments for help request (-h). If found,
    prints out help message. Needs to be passed a list of the options
    collected from the command line, as well as the help message.
    """

    import sys

    for a in options[:]:
        if a[0] == "-h":
            print message
            sys.exit(0)

