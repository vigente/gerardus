# Copyright (C) 2012 Dave Abrahams <dave@boostpro.com>
#
# Distributed under the Boost Software License, Version 1.0.  See
# accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt

import sys_path_setup
import glob, os, sys, itertools
from datetime import date, datetime
from warnings import warn
from subprocess import check_output
from xml.etree.cElementTree import ElementTree, Element
from dom import dashtag, xml_document, xmlns
from path import Path, pardir
import boost_metadata
from depgraph import *
from transitive import *
from uuid import uuid4 as make_uuid
from archive import Archive
from sign_feed import *
import threadpool
from read_dumps import read_dumps
_ = dashtag

os.environ['DISPLAY']=''

generated_feeds = set()

package_prefixes = ['Boost', 'Ryppl']
def split_package_prefix(package_name):
    for prefix in package_prefixes:
        n = len(prefix)
        if (len(package_name) > n
            and package_name.startswith(prefix)
            and package_name[n].isupper()):
            return prefix, package_name[n:]
    return None, package_name

def get_brand_name(cmake_name):
    prefix,base = split_package_prefix(cmake_name)
    return prefix + '.' + base if prefix else base

feed_uri_base = 'https://raw.github.com/ryppl/feeds/gh-pages/'
ryppl_feed_uri_base = feed_uri_base+'ryppl/'
boost_feed_uri_base = feed_uri_base+'boost/'

boost_icon = _.icon(
    href='http://svn.boost.org/svn/boost/website/public_html/live/gfx/boost-dark-trans.png', 
    type="image/png")

BSL_1_0 = 'OSI Approved :: Boost Software License 1.0 (BSL-1.0)'

def boost_feed_uri(name):
    return boost_feed_uri_base + name + '.xml'

def ryppl_feed_uri(name):
    return ryppl_feed_uri_base + name + '.xml'

class GenerateBoost(object):

    def cmake_package_to_feed_uri(self, cmake_package_name, component):
        prefix, basename = split_package_prefix(cmake_package_name)

        dump = self.dumps.get(cmake_package_name)
        repo = (dump.findtext('source-directory') - self.source_root) if dump else basename.lower()
        
        return feed_uri_base + '%s%s%s.xml' \
            % (
                  (prefix.lower() + '/' if prefix else '')
                , repo
                , ('' if component == 'bin' else '-dev' if prefix != 'Boost' else '-'+component)
              )

    def has_binary_lib(self, cmake_package_name):
        return self.dumps[cmake_package_name].find('libraries/library') is not None

    def _dev_requirements(self, cmake_package_names):
        return [
            _.requires(interface=self.cmake_package_to_feed_uri(cmake_package, 'dev')) [
                _.environment(insert='.', mode='replace', name=cmake_package+'_DIR')
                ]
            for cmake_package in set(cmake_package_names)
            ] +  [
            _.requires(interface=self.cmake_package_to_feed_uri(cmake_package, 'bin')) [
                # This one is needed for the library build
                _.environment(insert='.', mode='replace', name=cmake_package+'_BIN_DIR')
                ]
            for cmake_package in set(cmake_package_names) if self.has_binary_lib(cmake_package)
            ]

    def _run_requirements(self, cmake_package_names):
        return [
            _.requires(interface=self.cmake_package_to_feed_uri(cmake_package, 'bin')) [
                # These are needed for runtime (they come through compile:include-binary)
                _.environment(insert='bin', name='PATH')
              , _.environment(insert='lib', name='LD_LIBRARY_PATH')
              , _.environment(insert='lib', name='DYLD_LIBRARY_PATH')
                ]
            for cmake_package in set(cmake_package_names) if self.has_binary_lib(cmake_package)
            ]

    _empty_zipball = (
        _.archive(
            extract='empty', href=feed_uri_base+'empty.zip'
          , size=162, type="application/zip"
        ), 
        _.manifest_digest(sha1new='da39a3ee5e6b4b0d3255bfef95601890afd80709')
        )


    class GenerateRepo(object):
                
        def __getattr__(self, name):
            return getattr(self.ctx,name)

        preinstall_feed = None

        def __init__(self, ctx, cmake_name):
            self.ctx = ctx
            self.cmake_name = cmake_name
            self.in_cluster = self.cluster_map.get(cmake_name)
            dump = self.dumps[cmake_name]
            self.srcdir = dump.findtext('source-directory')
            self.git_revision = check_output(['git', 'rev-parse', 'HEAD'], cwd=self.srcdir).strip()
            self.repo = str(self.srcdir - self.source_root)
            self.boost_metadata = boost_metadata.lib_metadata(self.repo, self.boost_metadata)
            
            def dump_list(xpath):
                ret = [x.text for x in dump.findall(xpath)]
                return ret if len(ret)>0 else None

            self.binary_libs = dump_list('libraries/library')
            self.executables = dump_list('executables/executable')
            self.include_directories = dump_list('include-directories/directory')

            self.components = set()
            if self.binary_libs or self.executables:
                self.components.add('bin')
            if self.include_directories:
                self.components.add('dev')

            self.brand_name = get_brand_name(cmake_name)
            self.build_dependencies = self.transitive_dependencies[cmake_name]
            assert cmake_name not in self.build_dependencies, "Self-loop detected: "+cmake_name

            print '##', self.brand_name

            self.tasks.add_task(self._write_src_feed)

            if self.in_cluster:
                self.tasks.add_task(self._write_rawsrc_feed)
                self.preinstall_feed = self._cluster_feed_name(self.in_cluster)
                self.preinstall_subdirectory = self.repo + '/'

            # if there is more than one build product, we need a
            # preinstall feed.  Cluster preinstall feeds are handled
            # separately.
            elif len(self.components) > 1:
                self.preinstall_feed = self.repo + '-preinstall'
                self.preinstall_subdirectory = ''
                self.build('preinstall')

            if self.preinstall_feed:
                for c in self.components:
                    self.copy_from_preinstall(c)
            else:
                # No preinstall => header-only or executable-only, but not both
                sys.stdout.flush()
                assert len(self.components) <= 1
                for x in self.components:
                    self.build(iter(self.components).next())

        def copy_from_preinstall(self, component):
            self.tasks.add_task(
                self._write_feed, component
                , self._implementation('*-src') [ self._empty_zipball ]
                , _.command(name='compile') [
                    _.runner(interface=ryppl_feed_uri('ryppl')) [
                            _.arg['writable-directory-overlay'], _.arg['${SRCDIR}'], _.arg['${DISTDIR}']
                        ]
                    , _.requires(interface=boost_feed_uri(self.preinstall_feed)) [
                        _.environment(
                            insert=self.preinstall_subdirectory + component, 
                            mode='replace', name='SRCDIR')
                        ]
                    , self.implementation_template(component)
                    ]
                )
        
        def implementation_template(self, component):
            if component == 'dev' and not self.binary_libs:
                # Header-only -dev feeds can be marked architecture-independent
                # TODO: can we make -dev a *-* feed on Posix always, since there's no implib?
                return xmlns.compile.implementation(arch='*-*')
            elif component == 'bin' and self.executables:
                 # -bin feeds with executables should include their run commands
                return xmlns.compile.implementation() [
                          [ _.command(name='run', path='bin/'+x) for x in self.executables[:1] ]
                        , [ _.command(name=x, path='bin/'+x) for x in self.executables[1:]]
                          , self._run_requirements(self.build_dependencies)
                          ]
            else:
                return None
            
        def build(self, component):
            self.tasks.add_task(
                self._write_feed, component
                ,  self._git_snapshot('*-src') 
                     if component == 'preinstall' or not self.preinstall_feed 
                   else self._implementation('*-src')[ self._empty_zipball ]
                , _.command(name='compile') [
                    _.runner(interface=ryppl_feed_uri('ryppl')) [
                        _.arg[ '0install-cmake' ] 
                        , _.arg[ '--build-type=Debug' ] if component == 'dbg' else None
                        , [_.arg[ c ] for c in self.components]
                        ]

                    , self._dev_requirements(self.build_dependencies)
                    , self.implementation_template(component)
                    ]
                )

        
        def _feed_name(self, component):
            return self.repo + ('' if component == 'bin' else '-'+component) + '.xml'
            
        def _write_feed(self, component, *contents):
            interface = self._interface(component)[ 
                _.group(license=BSL_1_0) [
                    contents 
                ]
            ]
            interface.indent()
            feed_path = self.feed_dir/self._feed_name(component)
            assert feed_path not in generated_feeds
            generated_feeds.add(feed_path)
            xml_document(interface).write(feed_path, encoding='utf-8', xml_declaration=True)
            sys.stdout.flush()

            sign_feed(feed_path)
            

        def _feed_uri(self, component):
            return boost_feed_uri_base+self._feed_name(component)

        _human_component = {
            'bin':'binaries'
            , 'src':'source code'
            , 'rawsrc':'raw source code'
            , 'dev':'development files'
            , 'dbg':'debugging version'
            , 'preinstall':'built state'
            , 'doc':'documentation'
            }

        def _interface(self, component):
            interface = _.interface(
                uri=self._feed_uri(component)
              , xmlns='http://zero-install.sourceforge.net/2004/injector/interface'
              , **{
                    'xmlns:compile':'http://zero-install.sourceforge.net/2006/namespaces/0compile'
                  , 'xmlns:dc':'http://purl.org/dc/elements/1.1/'
                  , 'xmlns:ryppl':'http://ryppl.org/2012'
                    })[
                _.name['%s (%s)' % (self.brand_name, self._human_component[component])]
              , boost_icon
              , _('ryppl:vcs-repository')(
                    type="git", href='http://github.com/boost-lib/%s.git' % self.repo)
              ]

            # These tags can be dragged directly across from our lib_metadata
            for tag in 'summary','homepage','dc:author','description','category':
                interface <<= self.boost_metadata.findall(tag)

            return interface

        def _git_snapshot(self, arch):
            git_revision = check_output(['git', 'rev-parse', 'HEAD'], cwd=self.srcdir).strip()
            archive_uri = 'http://nodeload.github.com/boost-lib/' + self.repo + '/zipball/' + git_revision
            zipball = Archive(archive_uri, self.repo, git_revision)

            impl = self._implementation(arch) [
                _.archive(
                    extract=zipball.subdir
                    , href=archive_uri
                    , size=zipball.size
                    , type='application/zip')
              , _.manifest_digest(sha1new=zipball.digest)
            ]
            impl.element.attrib['ryppl:vcs-revision']=git_revision
            return impl

        def _write_src_feed(self):
            self_loop = set([self.cmake_name])

            dump = self.dumps[self.cmake_name]
            build_deps = (
                set(fp.findtext('arg') for fp in dump.findall('find-package'))
                | set(ud.text for ud in dump.findall('depends/dependency'))
                ) - self_loop

            src_dump = self.src_dumps[self.cmake_name]
            develop_deps = (
                set(fp.findtext('arg') for fp in src_dump.findall('find-package'))
                - build_deps - self_loop)

            self._write_feed(
                'src'
                , self._git_snapshot('*-*')
                , [
                    _.requires(interface=self.cmake_package_to_feed_uri(cmake_package, 'src'))
                    for cmake_package in build_deps
                    ]
                , _.command(name='develop') [
                    _.runner(interface=ryppl_feed_uri('ryppl'))
                    , [
                        _.requires(interface=self.cmake_package_to_feed_uri(cmake_package, 'src'))
                        for cmake_package in develop_deps
                        ]
                    ]
                )
                    
        def _write_rawsrc_feed(self):
            self._write_feed(
                'rawsrc'
               , self._git_snapshot('*-*')
                )

    def _build_dependencies(self, cmake_package_name):
        return self.transitive_dependencies[cmake_package_name]

    def _implementation(self, arch):
        return _.implementation(
            arch=arch, id=make_uuid(), released=date.today().isoformat(), 
            stability='developer', version=self.version)

    def _write_cluster_feed(self, cluster):
        feed_name = self._cluster_feed_name(cluster)
        
        subdirs = dict((self._src_dir_name(x),x) for x in cluster)
        
        interface = _.interface(
            uri=boost_feed_uri(feed_name)
            , xmlns='http://zero-install.sourceforge.net/2004/injector/interface'
            , **{
                'xmlns:compile':'http://zero-install.sourceforge.net/2006/namespaces/0compile'
                , 'xmlns:dc':'http://purl.org/dc/elements/1.1/'
                })[
            _.name['%s (built state)' % ', '.join(cluster),]
            , _.summary['Evil Build Dependency Cluster (EBDC)']
            , boost_icon
            , _.group(license=BSL_1_0)
            [
                self._implementation('*-src')[self._empty_zipball]
                , _.command(name='compile') [

                    _.runner(interface=ryppl_feed_uri('ryppl')) 
                    [
                        _.arg[ '0install-cmake' ]

                        # , _.arg[ '--build-type=Debug' ] if component == 'dbg' else []
                        
                        , [ _.arg[x] 
                            for d,c in subdirs.items()
                            for x in ['--add-subdirectory', '${%s_SRCDIR}' % c, d ]
                            ]

                        # Component selection may be too naive here in
                        # general, but hopefully this whole block will
                        # be obsolete as we eliminate clusters
                        , [ _.arg[ d+'/'+component ] for component in 'dev','bin' for d in subdirs ]
                        ]

                    , [
                        _.requires(interface=boost_feed_uri(d+'-rawsrc')) [
                            _.environment(insert='.', mode='replace', name=c+'_SRCDIR')
                            ]
                        for d,c in subdirs.items()]

                    , self._dev_requirements(
                            dep for x in cluster for dep in self._build_dependencies(x)
                            if dep not in cluster
                            )
                    ]
                ]
            ]

        interface.indent()
        feed_path = self.feed_dir/feed_name+'.xml'
        xml_document(interface).write(feed_path, encoding='utf-8', xml_declaration=True)
        sign_feed(feed_path)

    def _delete_old_feeds(self):
        print '### deleting old feeds...'
        for old_feed in glob.glob(os.path.join(self.feed_dir,'*.xml')):
            if Path(old_feed).name != 'CMakeLists.xml':
                os.unlink(old_feed)

    # I think this will end up being unused
    def _cluster_name(self, cluster):
        splits = [ split_package_prefix(x) for x in cluster ]
        prefix = splits[0][0]

        names = cluster
        if prefix and all(splits[1:][0] == prefix):
            names = [prefix] + [x[1] for x in splits]

        return '-'.join(names)
        
    def _src_dir_name(self, cmake_package_name):
        return Path(self.dumps[cmake_package_name].findtext('source-directory')).name

    def _cluster_feed_name(self, cluster):
        return '-'.join(self._src_dir_name(x) for x in cluster) + '-preinstall'

    def _find_dependency_cycles(self):
        print '### Checking for dependency cycles... ',
        from SCC import SCC

        successors = Successors(self.dumps)

        # Find all Strongly-Connected Components (SCCs) that contain
        # multiple vertices.  Each of these must be built as a unit
        self.clusters = set(
            tuple(sorted(s))
            for s in SCC(str,lambda v: successors(self.dumps,v)).getsccs(self.dumps)
            if len(s) > 1)
        
        if len(self.clusters) > 0:
            warn( 
                'Build dependency graph contains cycles.  All SCCs:\n'
                + repr(self.clusters))

        # Map each cmake module into its cluster
        self.cluster_map = dict(
            (lib, cluster) for cluster in self.clusters for lib in cluster)
        
    def __init__(self, dump_dir, feed_dir, source_root, site_metadata_file, src_dump_dir):
        self.dump_dir = dump_dir
        self.src_dump_dir = src_dump_dir
        self.feed_dir = feed_dir
        self.source_root = source_root

        self.dumps = read_dumps(self.dump_dir)
        self.src_dumps = read_dumps(self.src_dump_dir)

        self.transitive_dependencies = to_mutable_graph(self.dumps)
        inplace_transitive_closure(self.transitive_dependencies)
        
        # eliminate self-loops
        for v0,v1s in self.transitive_dependencies.items():
            v1s.discard(v0)
        
        # Make sure there are no modularity violations
        self._find_dependency_cycles()
        
        self.version = '1.49-post-' + datetime.utcnow().strftime("%Y.%m.%d.%H%M")
        print '### new version =', self.version

        print '### reading Boost library metadata...'
        t = ElementTree()
        t.parse(site_metadata_file)
        self.boost_metadata = t.getroot().findall('library')

        self._delete_old_feeds()
            
        use_threads = True
        class Tasks(object):
            if use_threads:
                def __init__(self):
                    self.pool = threadpool.ThreadPool(100)
                def add_task(self, f, *args):
                    #print '###', f.__name__, args[:1]
                    def tsk(*args):
                        #print '>>>', f.__name__, args[:1]
                        f(*args)
                        #print '<<<', f.__name__, args[:1]
                    self.pool.add_task(tsk, *args)
                def wait_completion(self):
                    self.pool.wait_completion()
            else:
                def add_task(self, f, *args): 
                    return f(*args)
                def wait_completion(self): 
                    pass

        self.tasks = Tasks()

        for cluster in self.clusters:
            self.tasks.add_task(self._write_cluster_feed, cluster)

        for cmake_name in self.dumps:
            self.GenerateRepo(self, cmake_name)

        print '### Awaiting completion...'
        self.tasks.wait_completion()
        print '### Done.'

if __name__ == '__main__':
    argv = sys.argv

    ryppl = Path('/Users/dave/src/ryppl')
    feeds = ryppl / 'feeds'
    lib_db_default = '/Users/dave/src/boost/svn/website/public_html/live/doc/libraries.xml'

    dump_dir = Path(argv[1]) if len(argv) > 1 else feeds/'dumps'
    GenerateBoost(
        dump_dir=dump_dir
      , feed_dir=Path(argv[2] if len(argv) > 2 else feeds/'boost')
      , source_root=Path(argv[3] if len(argv) > 3 else ryppl/'boost-zero'/'boost')
      , site_metadata_file=Path(argv[4] if len(argv) > 4 else lib_db_default)
      , src_dump_dir=Path(argv[5] if len(argv) > 5 else (dump_dir / pardir / 'src_dumps').normpath)
        )
