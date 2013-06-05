import sys_path_setup
import tempfile
import os
import shutil
from path import *

class TempDir(Path):
    def __new__(cls, delete=True, *args, **kw):
        ret = Path.__new__(cls, tempfile.mkdtemp(*args, **kw))
        ret.delete=delete
        ret.saved_wd = os.getcwd()
        return ret

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.__cleanup()

    def __cleanup(self):
        if self.delete:
            self.delete = False
            if os.getcwd() == self:
                os.chdir(self.saved_wd)
            try:
                shutil.rmtree(self)
            except:
                print 'failed to clean up', self
                pass # we don't care very much if it doesn't get cleaned up

    def __del__(self):
        self.__cleanup()



