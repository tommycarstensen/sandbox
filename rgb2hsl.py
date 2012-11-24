#!/software/bin/python
#Tommy Carstensen, University College Dublin, 2007

def main(R,G,B,):

    fR = float(R / 255.0)
    fG = float(G / 255.0)
    fB = float(B / 255.0)
    ma = max(fR, fG, fB)
    mi = min(fR, fG, fB)
    if ma == mi:
        H = 0
    elif ma == fR:
        H = (60 * ((fG - fB)/(ma - mi))) % 360
    elif ma == fG:
        H = 120 + 60 * ((fB - fR)/(ma - mi))
    elif ma == fB:
        H = 240 + 60 * ((fR - fG)/(ma - mi))
    L = (ma + mi) / 2
    if ma == mi:
        S = 0
    elif L <= 0.5:
        S = (ma - mi) / (2 * L)
    elif L > 0.5:
        S = (ma - mi) / (2 - 2 * L)

    return H, S, L

if __name__ == '__main__':
    main()
