def volumeDNA(volumefile):
    volume=open(volumefile).readlines()[0].strip('\n').split()
    DNA=float(volume[2])
    return DNA

def volumecomplex(volumefile):
    volume=open(volumefile).readlines()[0].strip('\n').split()
    DNA=float(volume[0])
    return DNA