import os, sys
sys.path = [
    os.path.join(
        os.path.dirname(
            os.path.dirname(
                os.path.abspath(__file__)))
        , 'ryppl'
        , 'support')
    ] + sys.path
