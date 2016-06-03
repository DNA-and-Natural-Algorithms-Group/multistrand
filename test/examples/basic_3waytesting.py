from multistrand.objects import *
from multistrand.options import Options
from multistrand.system import *

o = Options()
o.dangles = 2
o.temperature = 25



initialize_energy_model(o)

toe = Domain(name='toehold',sequence='TCTCCATGTC')
bm = Domain(name='branch',sequence='CCCTCATTCAATACCCTACG')
extra = Domain(name='extra',sequence='CCACATACATCATATT')

substrate = toe.C + bm.C
invasion = substrate.C
incumbent = extra + bm

c = Complex( strands=[incumbent,invasion,substrate],structure=".(+.(+))")
c = Complex( strands=[incumbent,invasion,substrate],structure = "................((((((((((((((((((((+....................((((((((((+))))))))))))))))))))))))))))))")

st = "................" + 20 * "(" + "+" + 20*"." + "((((((((((+))))))))))" + 20 * ")"
primary = [0] * 20
intermed = [0] * 19
for i in range(20):
    st = "................" + (20-i) * "(" + i*"." + "+" + (20-i)*"." + i*"(" + "((((((((((+))))))))))" + 20 * ")"
    c = Complex( strands=[incumbent,invasion,substrate],structure = st)
    primary[i] = energy([c],o)[0]
    
    if i>0:
        st_int = "................" + (20-i) * "(" + i*"." + "+" + (21-i)*"." + (i-1)*"(" + "((((((((((+))))))))))" + (i-1) * ")" + "." + (20-i) * ")"
        c = Complex( strands=[incumbent,invasion,substrate],structure = st_int)
        intermed[i-1] = energy([c],o)[0]

toeh = [0] * 6
for i in range(6):
    te = "................((((((((((((((((((((+....................((((" + (6-i)* "(" + i*"." + "+" + i*"." + (6-i)* ")" + "))))))))))))))))))))))))"
    c = Complex( strands=[incumbent,invasion,substrate],structure = te)
    toeh[i] = energy([c],o)[0]
    
toeh = toeh[::-1]
