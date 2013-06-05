# Substantially lifted from 0compile's autocompile.py
from zeroinstall.injector import arch, model, iface_cache
import platform
uname = arch._uname + platform.uname()[len(arch._uname):]

arch.machine_groups['missing'] = arch.machine_groups.get(uname[4], 0)
arch.machine_ranks['missing'] = max(arch.machine_ranks.values()) + 1
host_arch = '*-missing'

class MissingImplementation(model.ZeroInstallImplementation):
    # Assume that this (potential) binary is available so that we
    # can select it as a dependency.
    def is_available(self, stores):
        return True

class SloppyCache(iface_cache.IfaceCache):
    def __init__(self):
        iface_cache.IfaceCache.__init__(self)
        self.done = set()

    def get_feed(self, url, force = False):
        feed = iface_cache.IfaceCache.get_feed(self, url, force)
        if not feed: return None

        if feed not in self.done:
            self.done.add(feed)

            # For each implementation, add a corresponding "missing"
            # implementation

            srcs = [x for x in feed.implementations.itervalues() if x.arch]
            for x in srcs:
                new_id = 'ryppl=' + x.id
                if not new_id in feed.implementations:
                    new = MissingImplementation(feed, new_id, None)
                    feed.implementations[new_id] = new
                    new.set_arch(host_arch)
                    new.version = x.version

        return feed
