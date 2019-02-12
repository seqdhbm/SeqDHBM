class SeqResult():

    def __init__(self, seq, mode):
        self.wesa_mode = mode
        self.rawseq = seq
        self.line = seq
        self.header = ""
        self.result = {}
        self.analysis = """"""
        self.fail = False
        self.warnings = []


    def set_mode(self, m):
        self.wesa_mode = m
