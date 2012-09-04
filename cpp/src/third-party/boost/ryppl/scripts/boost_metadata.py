import re
from dom import xmlns, tag as _
from path import Path

def lib_metadata(source_subdirectory, all_metadata):
    """Build an ElementTree.Element representing the metadata we can
    extract from boost's libraries.xml for the library rooted in the
    given source subdirectory of Boost.

    In practice, source_subdirectory may be a parent directory of
    several sub-libraries (e.g. boost/math, with sub-libraries
    complex, special_functions, etc.).
"""

    # Gather all the libraries from Boost's library metadata
    libs = [
        l for l in all_metadata
        if l.findtext('key') == source_subdirectory
        or l.findtext('key').startswith(source_subdirectory+'/') ]
    
    # Find a library in the metadata whose <key> matches the source
    # subdirectory, if possible
    rootlib = ([ l for l in libs if l.findtext('key') == source_subdirectory ] + [None])[0]

    summary = 'A C++ library from Boost.org'
    homepage = 'http://www.boost.org/'
    if not rootlib is None:
        homepage += rootlib.findtext('documentation')

    # if there's more than one library, we have something for the <description> field.
    if len(libs) > 1:
        description = '\n\n'.join([l.findtext('key') + ':\n' + l.findtext('description') for l in libs])
    else:
        description = l.findtext('description')
    
    
    # Gather up all the distinct author names we can find
    authors = set()
    for l in libs:
        a0 = re.split('\s*(?:(?:,\s*(?:and\s+)?)|(?:and\s+))', l.findtext('authors'))
        for a in a0:
            authors.add(re.sub(r'\s+', ' ', a))

    # Gather up all the distinct categories
    categories = set()
    for l in libs:
        for c in l.findall('category'):
            categories.add(c.text)

    m = _.metadata[
        _.source_subdirectory[source_subdirectory],
        _.summary[summary],
        _.homepage[homepage],
        [xmlns.dc.author[a] for a in authors],
        [_.category(type='http://boost.org/lib-categories/2012')[c] for c in categories]
        ]

    if description:
        m <<= _.description[description]

    return m

if __name__ == '__main__':
    from xml.etree import cElementTree as ElementTree
    t = ElementTree.ElementTree()
    t.parse('/Users/dave/src/boost/svn/website/public_html/live/doc/libraries.xml')
    
    import glob
    libs_root = Path('/Users/dave/src/boost/svn/trunk/libs')
    for x in glob.glob( libs_root / '*' ):
        relpath = x - libs_root
        print 20*'#' + ' ' + relpath + ' ' + 20*'#'
        print lib_metadata(relpath, t.getroot().findall('library'))
    
