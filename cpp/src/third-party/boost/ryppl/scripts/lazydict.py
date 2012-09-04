class lazydict(dict):
    def __init__(self, factory, *args, **kw):
        self.__factory = factory
        dict.__init__(self, *args, **kw)

    __missing = object()
    def __getitem__(self, x):
        x = self.get(x,self.__missing)
        return self.__factory() if x is self.__missing else x

