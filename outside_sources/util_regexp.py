import re


class NamedRE_helper_tuples( dict ):
    unique_count = 0
    def __init__(self, *args, **kargs):
        if NamedRE_helper_tuples.unique_count > 0:
            raise ValueError("Only one helper tuples dict should be in existence.")
        
        NamedRE_helper_tuples.unique_count = 1
        super(NamedRE_helper_tuples,self).__init__()
        self.__tuple_dict__ = {}
        self.__init_tuple_flags__()


    def __init_tuple_flags__(self):
        self.__set_tuple_action__( "NamedGroup",
                                   lambda val, name: r"(?P<{name}>{value})".format( name=name, value=val),
                                   "A grouping which is assigned a specific name in the results.",
                                   no_direct = True,  # cannot be used in a sub-list tuple.
                                   group = True)
        

        self.__set_tuple_action__( "Group",
                                   lambda val,name=None: r"({value})".format(value=val),
                                   "A grouping without a specific name, but included in the results.",
                                   group=True )
        self.__set_tuple_action__( "LooseGroup",
                                    lambda val,name=None: r"(?:{value})".format(value=val),
                                   "A grouping for convenience, not included in results.",
                                   group=True
                                   )

        self.__set_tuple_action__( "OneOrMore",
                                   lambda val="",name=None: r"{value}+".format(value=val),
                                   "One or more copies of this element.",
                                   unary = True )

        self.__set_tuple_action__( "NoneOrMore",
                                   lambda val="",name=None: r"{value}*".format(value=val),
                                   "Any number of copies of this element, including none.",
                                   unary = True )
        self.__set_tuple_action__( "NoneOrOne",
                                   lambda val="",name=None: r"{value}?".format(value=val),
                                   "One copy of this element, or none.",
                                   unary = True )
        self.__set_tuple_action__( "NonGreedy",
                                   lambda val="",name=None: r"{value}?".format(value=val),
                                   "Modifies OneOrMore or NoneOrMore to be non-greedy.",
                                   unary = True,
                                   afterUnary = True)
        self.__set_tuple_action__( "Or",
                                   lambda val=None,name=None: r"|",
                                   "Either of the two adjacent expressions.",
                                   loose_or = True
                                   )
        
    def __set_tuple_action__(self,name,fn,doc,**kargs):
        # if len(args) < 3:
        #     raise ValueError("Must be called with at least 3 regular arguments: Name, Function and Docstring")
        self.__tuple_dict__[name] = fn    # function, usually a lambda, 2 named args.
        self.__tuple_dict__[name].name = name  # the function knows its own name.
        self.__tuple_dict__[name].doc = doc  # the function knows its docstring, too.
        self.__tuple_dict__[name].flags = {}   # while we're at it, it knows all the flags.
        for k,v in kargs.iteritems():
            self.__tuple_dict__[name].flags[k] = v
        # set all those flags as needed.

    def __getattr__(self,name):
        #        if name.startswith( ("  "))
        if self.__tuple_dict__.has_key(name):
            f = self.__tuple_dict__[name]
            return f
        else:
            raise AttributeError("No such attribute: {0}".format(name))


NamedRE_Dict = NamedRE_helper_tuples()
# Top level item 


class NamedRE_item( object ):
    """ Primary parent class of all NamedRE component items. Defines
    basic names and functions"""

    #So, why does this even exist?

    def __init__(self, *args, **kargs):
        super(NamedRE_item,self).__init__()
        # The owner field is necessary for doing any form of string
        # lookup.  This possibly makes the actual objects
        # non-portable, but that's always going to be an issue: since
        # strings must be interpreted in the context of the named RE
        # dictionary, if you only move some strings and not others,
        # you may change the overall language. The best idea would be
        # to use repr to get the objects so that they should
        # reconstruct everything automaticallyonce you build a new
        # NamedRE with them.
    # def __str__(self):
    #     super(NamedRE_item,self).__str__()
    # def __repr__(self):
    #     super(NamedRE_item,self).__repr__()
    def toBNF(self):
        pass

class NamedRE_string( object):
    def __init__(self, *args, **kargs ):
        if 'literal' in kargs:
            self.literal = kargs['literal']
            del kargs['literal']
        else:
            self.literal = False
        if 'owner' in kargs:
            self.owner = kargs['owner']
            del kargs['owner']
        else:
            self.owner = None

        self.data = str( args[0] )
        # Note: owner gets handled by the NamedRE_item handler before
        #  the str base class gets to process the args.

        # if this string is literal, never perform lookup via owner.
        
    def __str__(self):
        if self.literal:
            return self.data
        elif self.owner is None:
            return self.data
        elif isinstance(self.owner,dict):
            if self.data in self.owner:
                return str(self.owner[self.data])
            else:
                return self.data
        
    def __repr__(self):
        return r"NamedRE_string(r'{0}')".format( self.data )
    def toBNF(self):
        if self.literal:
            return r'"{0}"'.format( self.data )
        elif isinstance(self.owner,dict):
            if self.data in self.owner:
                return r'<{0}>'.format(self.data)
            else:
                return r'"{0}"'.format(self.data)
        else:
            return r'"{0}"'.format(self.data)
        

class NamedRE_list(list):
    def __init__(self, *args, **kargs ):
        if 'name' in kargs:
            self.name = kargs['name']
            del kargs['name']
        else:
            self.name = None
            
        super(NamedRE_list,self).__init__(*args, **kargs )
        # Note: owner gets handled by the NamedRE_item handler before
        #  the list base class gets to process the args.
        tmp_bools =  [v for v in self if isinstance(v,bool)]
        tmp_others = [i for i in self if not isinstance(i,bool)]


        self.grouping_flag = any( tmp_bools ) or \
                             all( tmp_bools ) and None or\
                             False

        self[:] = tmp_others

        if self.name is None and self.grouping_flag is True:
            raise ValueError("Must give a NamedRE_list a name.")

    def __str__(self):
        cval = r"".join(map(str,self))
        if self.grouping_flag == True:
            return r"(?P<{0}>{1})".format(self.name,cval)
        elif self.grouping_flag == False:
            return r"(?:{0})".format(cval)
        elif self.grouping_flag == None:
            return cval
        else:
            raise ValueError("Something went horribly wrong in NamedRE_list")
    def __repr__(self):
#        cval = map(repr,self)

        if self.grouping_flag is None:
            return r"NamedRE_list( {0}, name='{1}' )".format(cval, self.name )
        else:
            return r"NamedRE_list( {0}, name='{1}' )".format([self.grouping_flag] + self, self.name )

    def toBNF(self):
        pass


class NamedRE_tuple(list):
    def __init__(self, *args, **kargs ):
        if 'name' in kargs:
            self.name = kargs['name']
            del kargs['name']
        else:
            self.name = None
        if 'owner' in kargs:
            self.owner = kargs['owner']
            del kargs['owner']
        else:
            self.owner = None
            
        super(NamedRE_tuple,self).__init__( *args, **kargs )
        if len(self)>2:
            raise ValueError("NamedRE_tuple has at most 2 values.")
        

    def __str__(self):
        type_f = self[0]
        val =    self[1]

        if self.owner and val in self.owner:
            val = self.owner[val]
        
        if 'no_direct' in type_f.flags:
            if type_f.flags['no_direct'] is True and self.name is None:
                raise ValueError("Cannot use a named grouping tuple if it doesn't have a name!")
        
        if val is None:
            fl  =  'loose_or' in type_f.flags and type_f.flags['loose_or'] or \
                   'unary' in type_f.flags and type_f.flags['unary'] or \
                   False

            #TODO: no error checking for NonGreedy having to occur after a unary op yet.
            if fl is False:
                raise ValueError("Cannot use a tuple without a value, except if it's a loose or!")
            return type_f()
            #TODO: check for unary after non-group, possible warning?
        
        else:  # val has value!
            return type_f( str(val), self.name )


    def __repr__(self):
        return "NamedRE_tuple((NamedRE_Dict.{0},{1}), name={2})".format(
            self[0].name,
            repr(self[1]),
            repr(self.name)
            )

    def toBNF(self):
        pass

class NamedRE( dict ):
    """ NamedRE : base class dict

    This class helps generate regular expression strings with specific
    named components, essentially allowing you to translate a BNF into
    a regular expression for use in parsing a input line.

    Creation: can be created just like any other dictionary, though it expects all values to be either a string or a list of (string or bool).
    
    Boolean values indicate whether the corresponding key should be considered a 'named' value for regexp grouping (if True), a non-capturing group (if False), or not implicitly grouped at all (if no bool is in the list). Note that non-capturing groups are useful to guarantee composition of strings, but probably not a necessary feature.

    A string value, either in the value list or as the key's value itself, represents one of two things:
    1. An actual regular expression string component, such as r'\w+'.
    2. The name of another member of the dictionary, that should be looked up when evaluating the resulting expression for this value's 'key'.

    Example:
      names = NamedRE([
          ('first',[True,'name']),
          ('last',[True,'name']),
          ('name',r'[a-zA-Z]+'),
          ('title',[True,'dr',r'|','mr',r'|','ms',r'|','other']),
          ('dr',[False,r'Dr.',r'|',r'Doctor']),
          ('mr',[False,r'Mr.']),
          ('ms',[False,r'Mrs.',r'|',r'Miss',r'|',r'Ms.']),
          ('other',[False,r'Rev.',r'|',r'Fr.']),
          ('fullname',[True,'title',r'?\s*','first',r'?\s*','last']),
          ('fullname_only',[r'^\s*','fullname',r'\s*$'])
          ])
          
      fullname_re_string = names.get('fullname_only')
      print fullname_re_string
      
      --> '^\\s*(?P<fullname>(?P<title>(?:Dr.|Doctor)|(?:Mr.)|(?:Mrs.|Miss|Ms.)|(?:Rev.|Fr.))?\\s*(?P<first>[a-zA-Z]+)?\\s*(?P<last>[a-zA-Z]+))\\s*$'

      # Note that we used fullname_only so it has to match the entire
      # string. Could just use names.get('fullname') if we don't care
      # about surroundings.
      
      fullname_re = re.compile( fullname_re_string )
      nm = fullname_re.match('Mr. Joseph Schaeffer')
      print nm.groups()
      --> ('Mr. Joseph Schaeffer', 'Mr.', 'Joseph', 'Schaeffer')
      print nm.groupdict()
      --> {'first': 'Joseph',
           'fullname': 'Mr. Joseph Schaeffer',
           'last': 'Schaeffer',
           'title': 'Mr.'}
      nm2 = fullname_re.match('His Ridiculousness, Sir Zifnab the Humble')
      print nm2
      --> None
    """

    def __init__(self, *args, **kargs):
        super(NamedRE,self).__init__( *args, **kargs)
        for i in self.keys():
            self._load_key( i, None )

    def _load_key( self, key, initial ):
        val = None

        if isinstance(initial,bool):
            val = initial
        else:
            if initial is None and isinstance(self[key],list):
                internal_list = map( self._load_key,
                                     [None] * len(self[key]),
                                     self[key])
                val = NamedRE_list( internal_list,
                                    name=key)
            
            

            elif isinstance(initial, tuple) or key and isinstance(self[key],tuple):
                val = NamedRE_tuple( (self[key][0], self._load_key(None,self[key][1])) )
            elif isinstance(initial, str) or key and isinstance( self[key], str ):
                val = NamedRE_string( initial or self[key] )
            else:
                raise ValueError("Wrong value passed to load NamedRE")
            val.owner=self
            
        if key is None:
            return val
        else:
            self[key] = val
            
   
    def get_bnf( self, name ):
        pass


class NamedRE_old( dict ):
    """ NamedRE : base class dict

    This class helps generate regular expression strings with specific
    named components, essentially allowing you to translate a BNF into
    a regular expression for use in parsing a input line.

    Creation: can be created just like any other dictionary, though it expects all values to be either a string or a list of (string or bool).
    
    Boolean values indicate whether the corresponding key should be considered a 'named' value for regexp grouping (if True), a non-capturing group (if False), or not implicitly grouped at all (if no bool is in the list). Note that non-capturing groups are useful to guarantee composition of strings, but probably not a necessary feature.

    A string value, either in the value list or as the key's value itself, represents one of two things:
    1. An actual regular expression string component, such as r'\w+'.
    2. The name of another member of the dictionary, that should be looked up when evaluating the resulting expression for this value's 'key'.

    Example:
      names = NamedRE([
          ('first',[True,'name']),
          ('last',[True,'name']),
          ('name',r'[a-zA-Z]+'),
          ('title',[True,'dr',r'|','mr',r'|','ms',r'|','other']),
          ('dr',[False,r'Dr.',r'|',r'Doctor']),
          ('mr',[False,r'Mr.']),
          ('ms',[False,r'Mrs.',r'|',r'Miss',r'|',r'Ms.']),
          ('other',[False,r'Rev.',r'|',r'Fr.']),
          ('fullname',[True,'title',r'?\s*','first',r'?\s*','last']),
          ('fullname_only',[r'^\s*','fullname',r'\s*$'])
          ])
          
      fullname_re_string = names.get('fullname_only')
      print fullname_re_string
      
      --> '^\\s*(?P<fullname>(?P<title>(?:Dr.|Doctor)|(?:Mr.)|(?:Mrs.|Miss|Ms.)|(?:Rev.|Fr.))?\\s*(?P<first>[a-zA-Z]+)?\\s*(?P<last>[a-zA-Z]+))\\s*$'

      # Note that we used fullname_only so it has to match the entire
      # string. Could just use names.get('fullname') if we don't care
      # about surroundings.
      
      fullname_re = re.compile( fullname_re_string )
      nm = fullname_re.match('Mr. Joseph Schaeffer')
      print nm.groups()
      --> ('Mr. Joseph Schaeffer', 'Mr.', 'Joseph', 'Schaeffer')
      print nm.groupdict()
      --> {'first': 'Joseph',
           'fullname': 'Mr. Joseph Schaeffer',
           'last': 'Schaeffer',
           'title': 'Mr.'}
      nm2 = fullname_re.match('His Ridiculousness, Sir Zifnab the Humble')
      print nm2
      --> None
    """

    def __init__(self, *args, **kargs):
        super(NamedRE,self).__init__( *args, **kargs)
    
    def get( self, name ):
        if not self.has_key( name ):
            raise KeyError("No key {0} found in dictionary.".format(name))

        if isinstance(self[name],str):
            return self[name]

        if isinstance(self[name],tuple):
            return self._get_tuple( index_name=name, direct_tuple=None )
        if callable(self[name]):
            return self._get_tuple( index_name=name, direct_tuple=None )

        if isinstance(self[name],list):
            return self._get_list( name )
        
        raise ValueError("Unknown value type {0} for key {1}.".format( type(self[name]),name))

    def _get_list( self, name ):
        """ _get_list(self, name)
        Returns the regexp created from the list self[name]
         by iteration over each list member.
         Elements of type string are looked up directly, while
         tuples and lists are looked up by _get_tuple and _get_list
         recursively. Boolean values are a shorthand for some of the
         grouping functions which would otherwise be explicit flags
         for the tuple type.
         """

        grouping_flag = None
        cval =""
        for i in self[name]:
            if type(i) is bool:
                grouping_flag = i
            elif isinstance(i,str):
                cval = cval + self._dispatch(i)
            elif isinstance(i,tuple):
                cval = cval + self._get_tuple( index_name=None, direct_tuple=i )
            elif callable(i):
                cval = cval + self._get_tuple( index_name=None, direct_tuple=(i,))
            else:
                raise ValueError("Your lists have elements that are not strings, bools or tuples!")
                
        # once we have recursively processed all the members, compose
        # the resulting string with whatever grouping flags were set
        if grouping_flag is True:
            # named capturing group
            return r"(?P<{0}>{1})".format(name,cval)
        elif grouping_flag is False:
            # non-capturing group
            return r"(?:{0})".format(cval)
        else:
            return cval

    def _dispatch(self, i):
        # do a single step lookup on this element
        if i in self:
            # if it's in the dictionary and result is a string, replace with the contents
            if isinstance(self[i],str):
                return self[i]
            else:
                # contents were more complicated, call get recursively.
                return self.get(i)
        else:
            # The string isn't in the dictionary, so just use
            # the string as a literal.
            return i
        
   
    def get_bnf( self, name ):
        pass
        
                
class NamedRe_Test(object):
    """ Holds examples used to create regexps in reactions_obj.py and others. """
    
    # rxn_system = NamedRE( [("reactiontok",r"->"),
    #                   ("react_sep",  r"\+"),
    #                   ("prod_sep",   r"\+"),
    #                   ("react",      r"\w+"),
    #                   ("prod",       r"\w+"),
    #                   ("lhs", [True,"react",r"(?:\s*","react_sep",r"\s*","react",")*"]),
    #                   ("rhs", [True,"prod",r"(?:\s*","prod_sep",r"\s*","prod",")*"]),
    #                   ("reaction",[True,"lhs",r"?\s*","reactiontok",r"\s*","rhs",r"?"]),
    #                   ("ratesym","(\w+)"),
    #                   ("rate_1",[False,r"\(\s*Rate:?\s*","ratesym",r"\s*\)"]),
    #                   ("rate_2",[False,r"\(\s*","ratesym",r"\s*\)"]),
    #                   ("rate_3",[False,"ratesym"]),
    #                   ("rate",[True, "rate_1", r"|", "rate_2" , r"|", "rate_3"]),
    #                   ("fullreaction",[False,r"^\s*","reaction",r"\s*","rate",r"?\s*$"])
    #                   ])

    # rxn_system can express the following BNF, where
    # <reaction_expression> is rxn_system.get("fullreaction"), <symbol> is
    # one of rxn_system.get("ratesym"), rxn_system.get("prod"), or
    # rxn_system.get("react") as appropriate. The other pieces should
    # all be named obviously in the above dictionary.
    
    # <reaction_expression> ::= <lhs> "->" <rhs> <rate>
    # <lhs> ::= <symbol> "+" <lhs> | <symbol> | ""
    # <rhs> ::= <symbol> "+" <rhs> | <symbol> | ""
    # <rate> ::= "(" "Rate:" <symbol> ")" |
    #            "(" <symbol> ")" |
    #            <symbol>
    # <symbol> ::= [0-9a-zA-Z_]+
    #              ^^ here we abuse regexp notation. This is any
    #              sequence of length >= 1 composed entirely from
    #              characters within the []'s ('-' denotes ranges).
    # Note that "" denotes empty string, thus a LHS does not need to have any components!

    # Here's how we generate the whole string for use by re.compile:

    # reaction_expression = rxn_system.get("fullreaction")   # this is the one used in reactions_obj - Reaction class.
    # reaction_expression_re = re.compile( reaction_expression )

    # rate_expr = rxn_system.get("rate")
    # rate_expr_re = re.compile( rate_expr)


    # # the name example from the NamedRE doc string.
    # names = NamedRE([
    #     ('first',[True,'name']),
    #     ('last',[True,'name']),
    #     ('name',r'[a-zA-Z]+'),
    #     ('title',[True,'dr',r'|','mr',r'|','ms',r'|','other']),
    #     ('dr',[False,r'Dr.',r'|',r'Doctor']),
    #     ('mr',[False,r'Mr.']),
    #     ('ms',[False,r'Mrs.',r'|',r'Miss',r'|',r'Ms.']),
    #     ('other',[False,r'Rev.',r'|',r'Fr.']),
    #     ('fullname',[True,'title',r'?\s*','first',r'?\s*','last']),
    #     ('fullname_only',[r'^\s*','fullname',r'\s*$'])
    #     ])
      
    # fullname_re_string = names.get('fullname_only')
    # fullname_re = re.compile( fullname_re_string )


    
    

