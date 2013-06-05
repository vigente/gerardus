import os, tempfile, urllib2
import cPickle as pickle
from hashlib import md5
from subprocess import check_output

class Archive(object):
    cache_dir = '/tmp/ryppl_archive_cache'

    def __new__(cls, uri, source_subdir, revision):

        key = (uri,source_subdir,revision)
        cache_file = os.path.join(
            cls.cache_dir
          , md5(pickle.dumps(key)).hexdigest())

        try:
            arch = pickle.load(open(cache_file))
            if arch.key == key:
                return arch
        except Exception,e:
            pass

        result = object.__new__(cls, uri, source_subdir, revision)
        result.__init__(uri, source_subdir, revision)
        if not os.path.isdir(cls.cache_dir):
            os.makedirs(cls.cache_dir)
        pickle.dump(result, open(cache_file, 'w'))
        return result

    def __init__(self, uri, source_subdir, revision):
        self.key = (uri, source_subdir, revision)
        
        self.subdir = 'boost-lib-' + source_subdir + '-' + revision[:7]

        f = tempfile.NamedTemporaryFile(suffix='.zip')
        try:
            contents = urllib2.urlopen(uri).read()
        except:
            print 'failed uri:', uri
            import sys
            sys.stdout.flush()
            raise

        f.write(contents)
        f.flush()

        self.digest = check_output(
            ['0install', 'digest', '--algorithm=sha1new', f.name, self.subdir]
            ).strip().split('=')[1]

        self.size = len(contents)

if __name__ == '__main__':
    key = (
        'http://nodeload.github.com/boost-lib/timer/zipball/2a05e6539be22c036f16d069af7f40ae1e120c18'
      , 'timer'
      , '2a05e6539be22c036f16d069af7f40ae1e120c18')

    print 'first...'
    a = Archive(*key)
    print 'second...'
    b = Archive(*key)
    print 'done.'
