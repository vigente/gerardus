# Copyright David Abrahams 2007-2012. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
from xml.etree import cElementTree as ElementTree
from xml.etree.cElementTree import tostring

# This nasty hack is needed because cElementTree is imperfect.  Its
# "Element" is really just a function that returns an instance of an
# old-style class that isn't accessible by any other means.
_Element = type(ElementTree.Element(''))

from copy import deepcopy
from cStringIO import StringIO

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

class meta_xmlns(object.__class__):
    def __getattr__(self, name):
        return xmlns(name)

class xmlns(object):
    __metaclass__ = meta_xmlns
    def __init__(self, name):
        self.__name = name
    def __getattr__(self, tagname):
        return tag(self.__name + ':' + tagname)

class metatag(object.__class__):
    def __getattr__(self, name):
        return tag(name)

# For debugging purposes only
def dump(x, indent=0):
     if isinstance(x, tag):
         dump(x.element,indent)
     elif isinstance(x, _Element):
         print indent*' '+'_.'+x.tag, x.attrib, (x.text,x.tail)
         for y in x:
             dump(y, indent+2)
     else:
         print repr(x)

class tag(object):
    __metaclass__ = metatag

    def __init__(self, _tag_name, _extra={}, **_attributes):
        self.__element = ElementTree.Element(
            _tag_name
          , dict( [(k, unicode(v)) for k,v in _attributes.items()] )
          , **_extra
            )

    def __call__(self, **attributes):
        for k,v in attributes.items():
            if k.startswith('_'):
                k = k[1:]
            self.element.attrib[k] = unicode(v)
        return self

    def __ilshift__(self, x):
        self._flatten_append(self.element, x)
        return self
        
    def __getattr__(self, attr):
        return getattr(self.element, attr)

    @staticmethod
    def _flatten_append(l, x):
        if x is None:
            pass
        elif isinstance(x, tag):
            l.append(x.element)
        elif isinstance(x, _Element):
            l.append(x)
        elif isinstance(x,list) or isinstance(x,tuple):
            for y in x:
                tag._flatten_append(l, y)
        else:
            if len(l):
                l[-1].tail = (l[-1].tail or u'') + unicode(x)
            else:
                l.text =  (l.text or u'') + unicode(x)
            
            
    def __getitem__(self, new_children):
        tag._flatten_append(self.element, new_children)
        return self

    def __str__(self):
        x = deepcopy(self.element)
        indent(x)
        return tostring(x)

    def __repr__(self):
        return self.__class__.__name__ + ': ' + str(self)

    def __iter__(self):
        return iter(self._children)
    
    @property
    def element(self):
        return self.__element

    def indent(self):
        indent(self.element)

class dashmetatag(metatag):
    def __getattr__(self, name):
        return tag(name.replace('_','-'))
    
class dashtag(tag):
    __metaclass__ = dashmetatag
    

class PrettyElementTree(ElementTree.ElementTree):
    def __str__(self):
        x = deepcopy(self)
        indent(x.getroot())
        out = StringIO()
        x.write(out, encoding='utf-8', xml_declaration=True)
        return out.getvalue()

def xml_document(t):
    e = t.element if isinstance(t,tag) else t
    return PrettyElementTree(e)

if __name__ == '__main__':
    
    _ = tag
    print _.CarrierCode

    print _.RequestHeader[
                _.AccountNumber[ 33 ]
              , (_.MeterNumber[ 44 ]
              , _.CarrierCode[ 'FDXG' ])
            ]
    
    print xml_document(
            _.RequestHeader[
                _.AccountNumber[ 33 ]
              , _.MeterNumber[ 44 ]
              , _.CarrierCode[ 'FDXG' ]
            ])

    print xml_document(_('checkout-shopping-cart', xmlns="http://checkout.google.com/schema/2"))

    print xml_document(dashtag.checkout_shopping_cart(xmlns="http://checkout.google.com/schema/2"))
    
    t = _.text[ 'here is some ', 'text ', _.b['boldly ', 32], ' under rocks' ]
    print t
    assert str(t) == str(
        _.text[ None, 'here is some ', None, 'text ', _.b['boldly ', 32], ' under rocks', None ])

    r = _.text[ 'here is some ', 'text ', _.b['boldly ', 32], ' under rocks' ]
    r <<= _.i['fu']
    print r
    
    print xmlns.dc.name['somebody']
