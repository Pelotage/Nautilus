"""Granville similarity-law prototype (Python reference for crates/nautilus-ocean/src/fouling/granville.rs).

Original prototype from the research notes, unchanged except that the printed table now sits
under a main guard so derive_constants.py can import the functions. Run it directly to print
the friction-increase table the Rust test vectors were taken from.
"""
import math
K=0.42
nu=1.19e-6  # m2/s seawater 15C
def cf_s(Re):  # Schoenherr: 0.242/sqrt(CF)=log10(Re*CF)
    lo,hi=1e-4,0.02
    for _ in range(200):
        m=(lo+hi)/2
        f=0.242/math.sqrt(m)-math.log10(Re*m)
        if f>0: lo=m
        else: hi=m
    return m
def dU_nik(kp):
    if kp<3: return 0.0
    v=(1/K)*math.log(0.26*kp)
    if kp<15: v*=math.sin(math.pi/2*math.log10(kp/3)/math.log10(5))
    return max(v,0.0)
def dU_col(kp): return (1/K)*math.log(1+kp)
def cf_rough(Re_s, ks, L, dU):
    # outer: find L+ such that Re(L+) = Re_s
    def solve(Lp):
        kp=ks/L*Lp; d=dU(kp)
        # find CF where CF = cf_s(Re*exp(-K*d)) with Re = Lp/(sqrt(CF/2)*(1-sqrt(CF/2)/K))
        lo,hi=5e-4,0.03
        for _ in range(100):
            m=(lo+hi)/2
            u=math.sqrt(m/2); Re=Lp/(u*(1-u/K))
            g=m-cf_s(Re*math.exp(-K*d))
            if g>0: hi=m
            else: lo=m
        u=math.sqrt(m/2); return m, Lp/(u*(1-u/K))
    lo,hi=1e3,1e9
    for _ in range(80):
        Lp=math.sqrt(lo*hi); cf,Re=solve(Lp)
        if Re<Re_s: lo=Lp
        else: hi=Lp
    return cf
def ittc57(Re): return 0.075/(math.log10(Re)-2)**2
def ittc_dcf(ks,L,Re): return 0.044*((ks/L)**(1/3)-10*Re**(-1/3))+0.000125
if __name__ == "__main__":
    cases=[("124m frigate 15kn",124,15),("176m tanker 14kn",176,14),("230m KCS 24kn",230,24),("230m 16kn",230,16),("320m VLCC 15kn",320,15),("100m feeder 12kn",100,12)]
    ks_list=[30e-6,100e-6,150e-6,300e-6,1000e-6,3000e-6,10000e-6]
    for name,L,kn in cases:
        V=kn*0.5144; Re=V*L/nu; c0=cf_s(Re)
        print(f"\n{name}: Re={Re:.3e} CF_schoenherr={c0*1e3:.3f}e-3 CF_ITTC57={ittc57(Re)*1e3:.3f}e-3")
        for ks in ks_list:
            cn=cf_rough(Re,ks,L,dU_nik); cc=cf_rough(Re,ks,L,dU_col)
            it=ittc_dcf(ks,L,Re)
            print(f"  ks={ks*1e6:6.0f}um  dCF/CF Nikuradse-type={100*(cn-c0)/c0:6.1f}%  Colebrook-type={100*(cc-c0)/c0:6.1f}%  ITTC/Townsin dCF={it*1e3:.3f}e-3 ({100*it/ittc57(Re):5.1f}% of CF57)")
