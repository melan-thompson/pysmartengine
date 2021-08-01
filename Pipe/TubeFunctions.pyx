

cdef sgn(float x):
    if x>0:
        return 1
    elif x<0:
        return -1
    elif x==0:
        return 0

cdef minmod(float a, float b):
    if a*b<0:
        return 0
    else:
        sgn(a)*min(abs(a),abs(b))
